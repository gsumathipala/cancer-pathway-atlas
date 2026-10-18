(async function(){
"use strict";
let DATA;
try{ DATA=await (await fetch("data.json",{cache:"no-cache"})).json(); }
catch(e){ document.querySelector(".wrap").insertAdjacentHTML("afterbegin",'<p role="alert" style="padding:12px;border:1px solid var(--line);border-radius:8px">Could not load data.json. Serve this folder over HTTP (for example <code>python3 server.py</code>) rather than opening the file directly.</p>'); return; }
/* ---------- REFERENCES ---------- */
const REFS=DATA.refs;
const refLinks=a=>a.map(i=>`<a href="https://doi.org/${REFS[i-1][1]}" target="_blank" rel="noopener">[${i}]</a>`).join(" ");

/* ---------- PATHWAYS ---------- */
const TCGA10=DATA.order.tcga;
const OTHER=DATA.order.other;
const P=DATA.pathways;
const ALLP=[...TCGA10,...OTHER];
const pc=id=>`var(--p-${id})`;
function frStr(p){const f=p.fr; if(f.v!=null) return f.lab; if(f.lo!=null) return `≈${f.lo}–${f.hi}%`; return f.note||"";}
function frShort(p){const f=p.fr; if(f.v!=null) return f.v+"%"; if(f.lo!=null) return `${f.lo}–${f.hi}%`; return "";}

/* ---------- NODES ----------
 [id,label,pathway,role(O oncogene,T suppressor,D context-dependent,N rarely-mutated hub),freq%,x,y,alterations,note,cancers,drugs] */
const NODES=DATA.nodes.map(a=>({id:a[0],label:a[1],pw:a[2],role:a[3],f:a[4],x:a[5],y:a[6],alt:a[7],note:a[8],ca:a[9],dr:a[10]}));
const NOT_GENE=new Set(["AMPK","CASP","FUS"]);
const N=Object.fromEntries(NODES.map(n=>[n.id,n]));

/* ---------- EDGES  ">" activates, "|" inhibits, "-" complex ---------- */
const ESTR=DATA.edges;
const EDGES=ESTR.trim().split(/\s+/).map(s=>{const m=s.match(/^(\w+)([>|\-])(\w+)$/);return {s:m[1],t:m[3],k:m[2]===">"?"a":m[2]==="|"?"i":"b"}});
EDGES.forEach(e=>{if(!N[e.s]||!N[e.t]) console.warn("bad edge",e)});
const nbr={}; NODES.forEach(n=>nbr[n.id]=new Set());
EDGES.forEach(e=>{nbr[e.s].add(e.t);nbr[e.t].add(e.s)});

const ROLE={O:"Oncogene (gain of function)",T:"Tumour suppressor (loss of function)",D:"Context-dependent",N:"Pathway hub, rarely mutated"};
const rOf=f=>6.5+Math.sqrt(f)*2.6;

/* ---------- generic graph drawing ---------- */
function defs(svg,pfx){
  const d=svg.append("defs");
  [["a",""],["a","hi"],["i",""],["i","hi"]].forEach(([k,h])=>{
    const m=d.append("marker").attr("id",`${pfx}-${k}${h}`).attr("markerUnits","userSpaceOnUse").attr("orient","auto");
    const col=h?"var(--edge-hi)":"var(--edge)";
    if(k==="a"){m.attr("viewBox","0 0 10 10").attr("refX",9).attr("refY",5).attr("markerWidth",9).attr("markerHeight",9);
      m.append("path").attr("d","M0,0 L10,5 L0,10 z").style("fill",col);}
    else{m.attr("viewBox","0 0 4 14").attr("refX",2).attr("refY",7).attr("markerWidth",4).attr("markerHeight",13);
      m.append("rect").attr("width",4).attr("height",14).style("fill",col);}
  });
}
function edgePath(s,t,k){
  const dx=t.x-s.x,dy=t.y-s.y,L=Math.hypot(dx,dy)||1,ux=dx/L,uy=dy/L;
  const rs=s.r+2,rt=t.r+(k==="b"?2:4);
  return `M${s.x+ux*rs},${s.y+uy*rs} L${t.x-ux*rt},${t.y-uy*rt}`;
}
function drawGraph(g,pos,edges,pfx,opts){
  opts=opts||{};
  const eg=g.append("g").attr("class","edges");
  const es=eg.selectAll("path").data(edges).join("path")
    .attr("d",e=>edgePath(pos[e.s],pos[e.t],e.k))
    .style("fill","none").style("stroke","var(--edge)")
    .attr("stroke-width",e=>e.k==="b"?1.4:1.3)
    .attr("stroke-dasharray",e=>e.k==="b"?"2 3":(N[e.s].pw!==N[e.t].pw?"6 3":null))
    .attr("marker-end",e=>e.k==="b"?null:`url(#${pfx}-${e.k})`);
  const ng=g.append("g").attr("class","nodes");
  const ns=ng.selectAll("g").data(Object.values(pos)).join("g").attr("class","node")
    .attr("transform",d=>`translate(${d.x},${d.y})`).style("cursor","pointer")
    .attr("tabindex",0).attr("role","button").attr("aria-label",d=>`${N[d.id].label}, ${P[N[d.id].pw].short}`);
  ns.append("circle").attr("r",d=>d.r)
    .style("fill",d=>{const n=N[d.id];return n.role==="O"?pc(n.pw):"var(--surface)"})
    .style("stroke",d=>pc(N[d.id].pw))
    .attr("stroke-width",d=>{const r=N[d.id].role;return r==="T"?3.2:r==="D"?2.6:r==="N"?1.5:1.5})
    .attr("stroke-dasharray",d=>(N[d.id].role==="N"||d.ghost)?"3 2":null);
  ns.filter(d=>N[d.id].role==="D").append("path")
    .attr("d",d=>{const r=d.r-1.3;return `M0,${-r} A${r},${r} 0 0 1 0,${r} Z`}).style("fill",d=>pc(N[d.id].pw));
  ns.append("text").attr("y",d=>d.r+(opts.fs||12)+1).attr("text-anchor","middle")
    .style("font-family","var(--f-ui)").style("font-size",(opts.fs||12)+"px").style("fill","var(--fg)")
    .style("font-style",d=>NOT_GENE.has(d.id)?"normal":"italic").style("font-weight",500)
    .style("paint-order","stroke").style("stroke","var(--surface)").attr("stroke-width",3).attr("stroke-linejoin","round")
    .text(d=>N[d.id].label);
  if(opts.ghostOpacity) ns.filter(d=>d.ghost).style("opacity",opts.ghostOpacity);
  return {es,ns};
}

/* ---------- COMPOSITE MAP ---------- */
const svg=d3.select("#map");
defs(svg,"m");
const zoomG=svg.append("g");
const bg=zoomG.append("g");
bg.append("rect").attr("x",0).attr("y",0).attr("width",1700).attr("height",1090).style("fill","var(--surface)");
bg.append("rect").attr("x",14).attr("y",104).attr("width",1672).attr("height",972).attr("rx",36).style("fill","var(--cell)");
bg.append("rect").attr("x",14).attr("y",100).attr("width",1672).attr("height",16).attr("rx",8).style("fill","var(--mem)").style("opacity",.55);
bg.append("rect").attr("x",36).attr("y",655).attr("width",1628).attr("height",405).attr("rx",28).style("fill","var(--nuc)").style("stroke","var(--mem)").attr("stroke-width",4).attr("stroke-dasharray","14 6");
bg.append("rect").attr("x",1440).attr("y",280).attr("width",225).attr("height",280).attr("rx",60).style("fill","none").style("stroke","var(--mem)").attr("stroke-width",2.5).attr("stroke-dasharray","5 5");
const lab=(x,y,t,a)=>bg.append("text").attr("x",x).attr("y",y).attr("text-anchor",a||"start").style("font-family","var(--f-ui)").style("font-size","12px").style("letter-spacing",".14em").style("fill","var(--muted)").style("font-weight",600).text(t);
lab(24,30,"EXTRACELLULAR");lab(1676,92,"PLASMA MEMBRANE","end");lab(24,640,"CYTOPLASM");lab(52,685,"NUCLEUS");lab(1552,575,"MITOCHONDRION","middle");
// module captions
const cap=(x,y,t)=>bg.append("text").attr("x",x).attr("y",y).attr("text-anchor","middle").style("font-family","var(--f-display)").style("font-style","italic").style("font-size","13px").style("fill","var(--muted)").text(t);
cap(160,1062,"chromatin & epigenome");cap(380,1062,"splicing · telomere");cap(560,935,"MYC network");cap(820,1040,"RB checkpoint");cap(1110,1050,"p53 checkpoint");cap(1500,1040,"DNA repair");cap(1552,200,"apoptosis");
const posC={}; NODES.forEach(n=>posC[n.id]={id:n.id,x:n.x,y:n.y,r:rOf(n.f)});
const C=drawGraph(zoomG,posC,EDGES,"m");

const zoom=d3.zoom().scaleExtent([0.6,5]).translateExtent([[-200,-200],[1900,1290]])
  .filter(ev=>ev.type==="wheel"?(ev.ctrlKey||ev.metaKey):(!ev.button))
  .on("zoom",ev=>zoomG.attr("transform",ev.transform));
svg.call(zoom).on("dblclick.zoom",null);
d3.select("#zin").on("click",()=>svg.transition().duration(250).call(zoom.scaleBy,1.4));
d3.select("#zout").on("click",()=>svg.transition().duration(250).call(zoom.scaleBy,1/1.4));
d3.select("#zreset").on("click",()=>svg.transition().duration(300).call(zoom.transform,d3.zoomIdentity));

const tip=document.getElementById("tip"),mapbox=document.getElementById("mapbox");
C.ns.on("mouseenter",(ev,d)=>{const n=N[d.id];tip.textContent=`${n.label} · ${P[n.pw].short} · ≈${n.f}%`;tip.hidden=false;})
  .on("mousemove",ev=>{const b=mapbox.getBoundingClientRect();tip.style.left=(ev.clientX-b.left)+"px";tip.style.top=(ev.clientY-b.top)+"px";})
  .on("mouseleave",()=>{tip.hidden=true;})
  .on("click",(ev,d)=>{ev.stopPropagation();selectGene(d.id);})
  .on("keydown",(ev,d)=>{if(ev.key==="Enter"||ev.key===" "){ev.preventDefault();selectGene(d.id);}});
svg.on("click",()=>clearSel());

let focusPw=null,selGene=null;
function applyHighlight(){
  const litN=new Set(),litE=new Set();
  if(selGene){litN.add(selGene);nbr[selGene].forEach(x=>litN.add(x));EDGES.forEach((e,i)=>{if(e.s===selGene||e.t===selGene)litE.add(i)});}
  else if(focusPw){NODES.forEach(n=>{if(n.pw===focusPw)litN.add(n.id)});EDGES.forEach((e,i)=>{if(N[e.s].pw===focusPw||N[e.t].pw===focusPw){litE.add(i);litN.add(e.s);litN.add(e.t);}});}
  const any=selGene||focusPw;
  C.ns.style("opacity",d=>!any||litN.has(d.id)?1:.13);
  C.ns.select("circle").attr("stroke-width",function(d){const r=N[d.id].role;const base=r==="T"?3.2:r==="D"?2.6:1.5;return d.id===selGene?base+2:base});
  C.es.style("opacity",(e,i)=>!any||litE.has(i)?1:.07)
    .style("stroke",(e,i)=>any&&litE.has(i)?"var(--edge-hi)":"var(--edge)")
    .attr("marker-end",(e,i)=>e.k==="b"?null:`url(#m-${e.k}${any&&litE.has(i)?"hi":""})`);
  document.querySelectorAll("#chips .chip").forEach(c=>c.setAttribute("aria-pressed",String(c.dataset.pw===focusPw)));
}
function clearSel(){focusPw=null;selGene=null;applyHighlight();panelOverview();}
function selectGene(id,zoomTo){selGene=id;focusPw=null;applyHighlight();panelGene(id);
  if(zoomTo){const n=N[id],k=2.2;svg.transition().duration(500).call(zoom.transform,d3.zoomIdentity.translate(850-n.x*k,545-n.y*k).scale(k));}}
function focusPathway(pid){selGene=null;focusPw=focusPw===pid?null:pid;applyHighlight();focusPw?panelPathway(pid):panelOverview();}

// chips
const chips=d3.select("#chips");
ALLP.forEach(pid=>{const b=chips.append("button").attr("class","chip").attr("data-pw",pid).attr("aria-pressed","false").on("click",()=>focusPathway(pid));
  b.append("span").attr("class","sw").style("background",pc(pid));b.append("span").text(P[pid].short);});
// gene search
const dl=document.getElementById("genelist");
NODES.forEach(n=>{const o=document.createElement("option");o.value=n.label;dl.appendChild(o);});
const q=document.getElementById("q");
function findGene(v){v=v.trim().toUpperCase();if(!v)return null;
  return NODES.find(n=>n.label.toUpperCase()===v)||NODES.find(n=>n.id===v)||NODES.find(n=>n.label.toUpperCase().split(/[\s\/()–-]+/).includes(v))||NODES.find(n=>n.label.toUpperCase().includes(v));}
q.addEventListener("change",()=>{const n=findGene(q.value);if(n)selectGene(n.id,true);});
q.addEventListener("keydown",e=>{if(e.key==="Enter"){const n=findGene(q.value);if(n)selectGene(n.id,true);}});

// legend key
function legendKey(el){
  el.innerHTML=`
  <span><svg width="20" height="20"><circle cx="10" cy="10" r="7" style="fill:var(--fg)"/></svg> oncogene</span>
  <span><svg width="20" height="20"><circle cx="10" cy="10" r="6.5" style="fill:var(--surface);stroke:var(--fg)" stroke-width="3"/></svg> tumour suppressor</span>
  <span><svg width="20" height="20"><circle cx="10" cy="10" r="6.5" style="fill:var(--surface);stroke:var(--fg)" stroke-width="2.4"/><path d="M10,4 A6,6 0 0 1 10,16 Z" style="fill:var(--fg)"/></svg> context-dependent</span>
  <span><svg width="20" height="20"><circle cx="10" cy="10" r="6" style="fill:var(--surface);stroke:var(--fg)" stroke-width="1.5" stroke-dasharray="3 2"/></svg> hub, rarely mutated</span>
  <span><svg width="38" height="12"><line x1="2" y1="6" x2="28" y2="6" style="stroke:var(--fg)" stroke-width="1.4"/><path d="M27,1 L36,6 L27,11z" style="fill:var(--fg)"/></svg> activates</span>
  <span><svg width="38" height="14"><line x1="2" y1="7" x2="32" y2="7" style="stroke:var(--fg)" stroke-width="1.4"/><rect x="32" y="1" width="3.5" height="12" style="fill:var(--fg)"/></svg> inhibits</span>
  <span><svg width="38" height="12"><line x1="2" y1="6" x2="36" y2="6" style="stroke:var(--fg)" stroke-width="1.4" stroke-dasharray="2 3"/></svg> complex</span>
  <span><svg width="38" height="12"><line x1="2" y1="6" x2="36" y2="6" style="stroke:var(--fg)" stroke-width="1.4" stroke-dasharray="6 3"/></svg> cross-pathway</span>
  <span>node area ∝ pan-cancer alteration frequency</span>`;
}
legendKey(document.getElementById("legendkey"));legendKey(document.getElementById("legendkey2"));

// panel
const panel=document.getElementById("panel");
function panelOverview(){
  const top=[...TCGA10].slice(0,5);
  panel.innerHTML=`<div class="eyebrow">How to read the map</div>
  <h2>One cell, every common lesion</h2>
  <p class="lead">Each node is a driver gene placed where its protein works: membrane receptors on top, cytoplasmic transducers in the middle, transcription, cell-cycle, repair and chromatin machinery in the nucleus. Colour gives the pathway, fill gives how the gene is broken, and size gives how often.</p>
  <p class="lead">Two signalling axes, RTK–RAS–MAPK and PI3K–mTOR, converge on MYC and cyclin D. Two checkpoints, RB and p53, are the brakes that almost every carcinoma disables, often both at once through <span class="gene">CDKN2A</span>.</p>
  <div><div class="eyebrow" style="margin-bottom:6px">Most often altered pathways</div>
  ${top.map(p=>`<div style="display:flex;justify-content:space-between;gap:8px;border-bottom:1px solid var(--line);padding:5px 0"><button class="chip" data-open="${p}"><span class="sw" style="background:${pc(p)}"></span>${P[p].name}</button><span class="freqnum">${frShort(P[p])}</span></div>`).join("")}</div>
  <p class="refs-inline">89% of 9,125 TCGA tumours carry at least one driver alteration in the ten canonical pathways ${refLinks([1])}.</p>`;
  panel.querySelectorAll("[data-open]").forEach(b=>b.addEventListener("click",()=>focusPathway(b.dataset.open)));
}
function panelPathway(pid){
  const p=P[pid],genes=NODES.filter(n=>n.pw===pid).sort((a,b)=>b.f-a.f);
  panel.innerHTML=`<div><span class="pill" style="color:${pc(pid)}">${TCGA10.includes(pid)?"TCGA canonical pathway":"Additional pathway"}</span></div>
  <h2>${p.name}</h2>
  <div class="freqnum" style="font-size:13px;color:var(--muted)">${frStr(p)}</div>
  <p class="lead">${p.mech}</p>
  <div><div class="eyebrow" style="margin-bottom:6px">Genes in this pathway</div><div class="nbrs">${genes.map(g=>`<button data-g="${g.id}">${g.label}</button>`).join("")}</div></div>
  <button class="btn" id="openpw">Open pathway view</button>
  <p class="refs-inline">${refLinks(p.refs)}</p>`;
  panel.querySelectorAll("[data-g]").forEach(b=>b.addEventListener("click",()=>selectGene(b.dataset.g)));
  panel.querySelector("#openpw").addEventListener("click",()=>{showPathway(pid);setView("pathways");});
}
function panelGene(id){
  const n=N[id],p=P[n.pw];
  const nb=[...nbr[id]].map(x=>N[x]);
  panel.innerHTML=`<div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap"><span class="pill" style="color:${pc(n.pw)}">${p.short}</span><span class="pill" style="color:var(--muted)">${ROLE[n.role]}</span></div>
  <h2 class="${NOT_GENE.has(id)?"":"gene"}">${n.label}</h2>
  <p class="lead">${n.note.replace(/\[(\d+(?:,\s*\d+)*)\]/g,(m,g)=>refLinks(g.split(/,\s*/).map(Number)))}</p>
  <dl class="kv">
    <dt>Frequency</dt><dd class="freqnum">≈${n.f}% of tumours, pan-cancer</dd>
    ${n.alt?`<dt>Lesions</dt><dd>${n.alt}</dd>`:""}
    ${n.ca?`<dt>Enriched in</dt><dd>${n.ca}</dd>`:""}
    ${n.dr?`<dt>Drugs</dt><dd>${n.dr}</dd>`:""}
  </dl>
  ${nb.length?`<div><div class="eyebrow" style="margin-bottom:6px">Direct partners</div><div class="nbrs">${nb.map(g=>`<button data-g="${g.id}" style="border-color:${pc(g.pw)}">${g.label}</button>`).join("")}</div></div>`:""}
  <button class="btn" id="openpw">Open ${p.short} pathway</button>`;
  panel.querySelectorAll("[data-g]").forEach(b=>b.addEventListener("click",()=>selectGene(b.dataset.g)));
  panel.querySelector("#openpw").addEventListener("click",()=>{showPathway(n.pw);setView("pathways");});
}
panelOverview();

/* ---------- INDIVIDUAL PATHWAYS ---------- */
const pwsvg=d3.select("#pwsvg");
const pwlist=d3.select("#pwlist");
let curPw="rtk";
[["Ten canonical (TCGA)",TCGA10],["Additional circuits",OTHER]].forEach(([g,ids])=>{
  pwlist.append("div").attr("class","eyebrow grp").text(g);
  ids.forEach(pid=>{const b=pwlist.append("button").attr("role","tab").attr("data-pw",pid).attr("aria-selected","false").on("click",()=>showPathway(pid));
    b.append("span").attr("class","sw").style("background",pc(pid));b.append("span").text(P[pid].short);b.append("span").attr("class","fq").text(frShort(P[pid]));});
});
function showPathway(pid){
  curPw=pid;
  pwlist.selectAll("button").attr("aria-selected",function(){return String(this.dataset.pw===pid)});
  const members=NODES.filter(n=>n.pw===pid).map(n=>n.id);
  const mset=new Set(members),ghosts=new Set();
  const es=EDGES.filter(e=>mset.has(e.s)||mset.has(e.t));
  es.forEach(e=>{if(!mset.has(e.s))ghosts.add(e.s);if(!mset.has(e.t))ghosts.add(e.t);});
  const all=[...members,...ghosts].map(id=>N[id]).sort((a,b)=>a.y-b.y||a.x-b.x);
  const rows=[];let cur=null;
  all.forEach(n=>{if(!cur||n.y-cur.y0>45){cur={y0:n.y,items:[]};rows.push(cur);}cur.items.push(n);});
  rows.forEach(r=>r.items.sort((a,b)=>a.x-b.x));
  const W=900,rowH=112,top=56,H=top+rows.length*rowH+20;
  pwsvg.attr("viewBox",`0 0 ${W} ${H}`).attr("aria-label",`${P[pid].name} pathway diagram`);
  pwsvg.selectAll("*").remove();
  defs(pwsvg,"pw");
  pwsvg.append("rect").attr("width",W).attr("height",H).style("fill","var(--surface)");
  pwsvg.append("text").attr("x",20).attr("y",30).style("font-family","var(--f-display)").style("font-size","20px").style("font-weight",600).style("fill",pc(pid)).text(P[pid].name);
  pwsvg.append("text").attr("x",W-20).attr("y",30).attr("text-anchor","end").style("font-family","var(--f-mono)").style("font-size","12px").style("fill","var(--muted)").text(frStr(P[pid]));
  const pos={};
  rows.forEach((r,ri)=>{const n=r.items.length;const y=top+ri*rowH+rowH/2-10;
    const comp=r.y0<130?"membrane":r.y0<650?"cytoplasm":"nucleus";
    if(ri===0||((rows[ri-1].y0<130?"membrane":rows[ri-1].y0<650?"cytoplasm":"nucleus")!==comp)){
      pwsvg.append("text").attr("x",20).attr("y",y-rowH/2+18).style("font-family","var(--f-ui)").style("font-size","10.5px").style("letter-spacing",".14em").style("font-weight",600).style("fill","var(--muted)").text(comp.toUpperCase());
      if(ri>0) pwsvg.append("line").attr("x1",20).attr("x2",W-20).attr("y1",y-rowH/2+4).attr("y2",y-rowH/2+4).style("stroke","var(--line)").attr("stroke-dasharray","4 4");
    }
    r.items.forEach((nd,i)=>{const ghost=!mset.has(nd.id);pos[nd.id]={id:nd.id,x:60+(W-120)*(n===1?.5:i/(n-1)),y:y+(n>4&&i%2?14:0),r:ghost?Math.min(rOf(nd.f),12):Math.max(rOf(nd.f)*1.25,12),ghost};});
  });
  const G=drawGraph(pwsvg.append("g"),pos,es,"pw",{fs:13,ghostOpacity:.45});
  G.ns.on("click",(ev,d)=>{setView("map");selectGene(d.id,true);})
    .on("keydown",(ev,d)=>{if(ev.key==="Enter"){setView("map");selectGene(d.id,true);}});
  G.ns.append("title").text(d=>`${N[d.id].label}: ≈${N[d.id].f}%`);
  const p=P[pid],genes=members.map(id=>N[id]).sort((a,b)=>b.f-a.f);
  document.getElementById("pwtext").innerHTML=`
    <div><div class="eyebrow">${TCGA10.includes(pid)?"TCGA canonical pathway":"Additional circuit"}</div><h2>${p.name}</h2>
    <div class="freqnum" style="color:var(--muted);font-size:13px">${frStr(p)}${p.fr.src?" "+refLinks([p.fr.src]):p.fr.ver?" "+refLinks([1]):""}</div></div>
    <section><div class="eyebrow">Normal circuit, in brief</div><p>${p.mech}</p></section>
    ${tutorialHTML(p)}
    <section><div class="eyebrow">How cancers break it</div><p>${linkify(p.lesion)}</p></section>
    <section><div class="eyebrow">Where it dominates</div><p>${linkify(p.where)}</p></section>
    <section><div class="eyebrow">Vulnerabilities</div><p>${linkify(p.vuln)}</p></section>
    ${evidenceHTML(p)}
    ${therapyHTML(p)}
    <section><div class="eyebrow">Most altered members</div><div class="nbrs">${genes.map(g=>`<button data-g="${g.id}">${g.label} <span class="freqnum" style="font-style:normal;color:var(--muted)">${g.f}%</span></button>`).join("")}</div></section>
    <p class="refs-inline">Sources ${refLinks(p.refs)}</p>`;
  document.querySelectorAll("#pwtext [data-g]").forEach(b=>b.addEventListener("click",()=>{setView("map");selectGene(b.dataset.g,true);}));
}
function linkify(s){return s.replace(/\[(\d+(?:,\s*\d+)*)\]/g,(m,g)=>refLinks(g.split(/,\s*/).map(Number)));}
const esc=s=>String(s).replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
function statusPill(t){if(t[1]==="—")return '<span class="st st-n">None</span>';
  return t[3]==="A"?`<span class="st st-a">FDA${t[4]?" "+t[4]:""}</span>`:t[3]==="T"?'<span class="st st-t">Trials</span>':'<span class="st st-x">Stopped</span>';}
function tutorialHTML(p){if(!p.steps)return "";
  return `<section><div class="eyebrow">How it works, step by step</div><ol class="steps">${p.steps.map((s,i)=>`<li><span class="sn">${i+1}</span><div><h4>${s[0]}</h4><p>${s[1]}${s[2]&&s[2].length?" "+refLinks(s[2]):""}</p></div></li>`).join("")}</ol></section>`;}
function evidenceHTML(p){if(!p.evidence)return "";
  return `<section><div class="eyebrow">Key evidence</div><ol class="timeline">${p.evidence.map(e=>`<li><span class="yr">${e[0]}</span><p>${e[1]} ${refLinks(e[2])}</p></li>`).join("")}</ol></section>`;}
function therapyHTML(p){if(!p.therapies)return "";
  const rows=p.therapies;
  return `<section><div class="eyebrow">Targeted therapies</div><div class="scroll"><table class="ther"><thead><tr><th>Agent</th><th>Target</th><th>Biomarker / setting</th><th>Status</th></tr></thead><tbody>${rows.map(t=>t[1]==="—"?`<tr><td colspan="3"><em>${esc(t[0])}</em></td><td>${statusPill(t)}</td></tr>`:`<tr><td><b>${esc(t[0])}</b>${t[5]&&t[5].length?" "+refLinks(t[5]):""}</td><td>${esc(t[1])}</td><td>${esc(t[2])}</td><td>${statusPill(t)}</td></tr>`).join("")}</tbody></table></div><p class="hint">FDA year is the first US approval of that agent in this pathway context; labels and indications change, so check the current prescribing information.</p></section>`;}
showPathway("rtk");

/* ---------- FREQUENCY ---------- */
const STATS=DATA.stats;
document.getElementById("stats").innerHTML=STATS.map(([b,t,r])=>`<div class="stat"><b>${b}</b><span>${t} ${refLinks([r])}</span></div>`).join("");
(function(){
  const s=d3.select("#chart-pw"),W=640,H=380,m={l:118,r:110,t:14,b:34};
  const x=d3.scaleLinear().domain([0,60]).range([m.l,W-m.r]);
  const order=[...TCGA10].sort((a,b)=>{const f=p=>P[p].fr.v!=null?P[p].fr.v:(P[p].fr.lo+P[p].fr.hi)/2;return f(b)-f(a)});
  const y=d3.scaleBand().domain(order).range([m.t,H-m.b]).padding(.32);
  [0,10,20,30,40,50,60].forEach(v=>{s.append("line").attr("x1",x(v)).attr("x2",x(v)).attr("y1",m.t).attr("y2",H-m.b).style("stroke","var(--line)");
    s.append("text").attr("x",x(v)).attr("y",H-m.b+16).attr("text-anchor","middle").style("font","11px var(--f-mono)").style("fill","var(--muted)").text(v+"%");});
  order.forEach(pid=>{const f=P[pid].fr,yy=y(pid),bh=y.bandwidth();
    s.append("text").attr("x",m.l-10).attr("y",yy+bh/2+4).attr("text-anchor","end").style("font","500 13px var(--f-ui)").style("fill","var(--fg)").text(P[pid].short);
    if(f.v!=null){s.append("rect").attr("x",x(0)).attr("y",yy).attr("width",x(f.v)-x(0)).attr("height",bh).attr("rx",2).style("fill",pc(pid));
      s.append("text").attr("x",x(f.v)+6).attr("y",yy+bh/2+4).style("font","500 12px var(--f-mono)").style("fill","var(--fg)").text(f.v+"% (stated)");}
    else{s.append("rect").attr("x",x(0)).attr("y",yy+bh*.3).attr("width",x(f.lo)-x(0)).attr("height",bh*.4).style("fill",pc(pid)).style("opacity",.35);
      s.append("rect").attr("x",x(f.lo)).attr("y",yy).attr("width",x(f.hi)-x(f.lo)).attr("height",bh).attr("rx",2).style("fill",pc(pid)).style("opacity",.75);
      s.append("text").attr("x",x(f.hi)+6).attr("y",yy+bh/2+4).style("font","12px var(--f-mono)").style("fill","var(--muted)").text(`≈${f.lo}–${f.hi}%`);}
  });
  s.append("text").attr("x",x(30)).attr("y",H-4).attr("text-anchor","middle").style("font","11px var(--f-ui)").style("fill","var(--muted)").text("share of tumours with ≥1 driver alteration in the pathway");
})();
(function(){
  const genes=[...NODES].filter(n=>n.role!=="N").sort((a,b)=>b.f-a.f).slice(0,25);
  const s=d3.select("#chart-genes"),W=640,H=640,m={l:150,r:60,t:10,b:34};
  const x=d3.scaleLinear().domain([0,40]).range([m.l,W-m.r]);
  const y=d3.scaleBand().domain(genes.map(g=>g.id)).range([m.t,H-m.b]).padding(.28);
  [0,10,20,30,40].forEach(v=>{s.append("line").attr("x1",x(v)).attr("x2",x(v)).attr("y1",m.t).attr("y2",H-m.b).style("stroke","var(--line)");
    s.append("text").attr("x",x(v)).attr("y",H-m.b+16).attr("text-anchor","middle").style("font","11px var(--f-mono)").style("fill","var(--muted)").text(v+"%");});
  genes.forEach(g=>{const yy=y(g.id),bh=y.bandwidth();
    const row=s.append("g").style("cursor","pointer").on("click",()=>{setView("map");selectGene(g.id,true);});
    row.append("text").attr("x",m.l-10).attr("y",yy+bh/2+4).attr("text-anchor","end").style("font","italic 500 12.5px var(--f-ui)").style("fill","var(--fg)").text(g.label);
    row.append("rect").attr("x",x(0)).attr("y",yy).attr("width",Math.max(2,x(g.f)-x(0))).attr("height",bh).attr("rx",2).style("fill",pc(g.pw));
    row.append("text").attr("x",x(g.f)+6).attr("y",yy+bh/2+4).style("font","11.5px var(--f-mono)").style("fill","var(--muted)").text("≈"+g.f+"%");
  });
  s.append("text").attr("x",x(20)).attr("y",H-4).attr("text-anchor","middle").style("font","11px var(--f-ui)").style("fill","var(--muted)").text("approximate share of all tumours, pan-cancer");
})();
(function(){
  const cols=["rtk","pi3k","cc","p53","myc","wnt","notch","hippo","tgfb","nrf2"];
  const CAN=DATA.cancers;
  const tl=["<5%","5–15%","15–30%","30–60%",">60%"];
  document.getElementById("heatkey").innerHTML=tl.map((t,i)=>`<span><i class="t${i}"></i>${t}</span>`).join("")+"<span>of tumours altered</span>";
  const tb=document.getElementById("heat");
  tb.innerHTML=`<thead><tr><th></th>${cols.map(c=>`<th><span class="sw" style="display:inline-block;width:8px;height:8px;border-radius:50%;background:${pc(c)};margin-right:4px"></span>${P[c].short}</th>`).join("")}</tr></thead>
  <tbody>${CAN.map(([code,name,t])=>`<tr><th class="row">${name}<small>${code}</small></th>${t.map((v,i)=>`<td class="t${v}" title="${name}: ${P[cols[i]].short} ${tl[v]}">${tl[v]}</td>`).join("")}</tr>`).join("")}</tbody>`;
})();

/* ---------- HALLMARKS ---------- */
(function(){
  const H=DATA.hallmarks;
  const M=DATA.hallmarkMap;
  const t=document.getElementById("hm");
  t.innerHTML=`<thead><tr><th></th>${H.map(([k,l],i)=>`<th class="${i>=10?"grpE":""}"><div>${l}</div></th>`).join("")}</tr></thead><tbody>${ALLP.map(pid=>{
    const m=Object.fromEntries(M[pid].split(" ").map(s=>s.split(":")));
    return `<tr><th class="row"><span class="sw" style="display:inline-block;width:10px;height:10px;border-radius:50%;background:${pc(pid)};margin-right:6px"></span>${P[pid].name}</th>${H.map(([k,l])=>{const v=m[k];
      return `<td title="${P[pid].short}: ${l}${v?v==="2"?" (primary)":" (secondary)":""}">${v==="2"?`<span class="dot" style="background:${pc(pid)}"></span>`:v==="1"?`<span class="dot" style="border:2px solid ${pc(pid)}"></span>`:""}</td>`}).join("")}</tr>`}).join("")}</tbody>`;
})();

/* ---------- PRINCIPLES ---------- */
const PR=DATA.principles;
document.getElementById("principles").innerHTML=PR.map(([h,t],i)=>`<li><span class="n">${String(i+1).padStart(2,"0")}</span><div><h3>${h}</h3><p>${linkify(t)}</p></div></li>`).join("");
document.getElementById("refs").innerHTML=REFS.map(([c,d])=>`<li>${c}. <a class="doi" href="https://doi.org/${d}" target="_blank" rel="noopener">doi:${d}</a></li>`).join("");

/* ---------- THERAPIES ---------- */
(function(){
  const rows=[];ALLP.forEach(pid=>(P[pid].therapies||[]).forEach(t=>{if(t[1]!=="—")rows.push({pid,t});}));
  const sel=document.getElementById("th-pw"),qq=document.getElementById("th-q");
  ALLP.forEach(pid=>{const o=document.createElement("option");o.value=pid;o.textContent=P[pid].name;sel.appendChild(o);});
  let st="all";
  document.querySelectorAll("#th-status button").forEach(b=>b.addEventListener("click",()=>{st=b.dataset.st;document.querySelectorAll("#th-status button").forEach(x=>x.setAttribute("aria-pressed",String(x===b)));render();}));
  sel.addEventListener("change",render);qq.addEventListener("input",render);
  function render(){
    const pv=sel.value,qv=qq.value.trim().toLowerCase();
    const f=rows.filter(r=>(pv==="all"||r.pid===pv)&&(st==="all"||r.t[3]===st)&&(!qv||r.t.slice(0,3).join(" ").toLowerCase().includes(qv)));
    const appr=rows.filter(r=>r.t[3]==="A").length;
    document.getElementById("th-count").textContent=`${f.length} of ${rows.length} entries · ${appr} approved agent groups across ${new Set(rows.filter(r=>r.t[3]==="A").map(r=>r.pid)).size} pathways`;
    document.getElementById("th-body").innerHTML=f.map(({pid,t})=>`<tr><td><button class="chip" data-open="${pid}"><span class="sw" style="background:${pc(pid)}"></span>${P[pid].short}</button></td><td><b>${esc(t[0])}</b></td><td>${esc(t[1])}</td><td>${esc(t[2])}</td><td>${statusPill(t)}</td><td class="refs-inline">${t[5]&&t[5].length?refLinks(t[5]):""}</td></tr>`).join("")||'<tr><td colspan="6" class="hint">No therapies match these filters.</td></tr>';
    document.querySelectorAll("#th-body [data-open]").forEach(b=>b.addEventListener("click",()=>{showPathway(b.dataset.open);setView("pathways");}));
  }
  render();
})();

/* ---------- VIEWS ---------- */
const VIEWS=["map","pathways","frequency","therapies","hallmarks","principles"];
function setView(v){if(!VIEWS.includes(v))v="map";
  VIEWS.forEach(x=>{document.getElementById("view-"+x).hidden=x!==v;document.getElementById("tab-"+x).setAttribute("aria-selected",String(x===v));});
  try{history.replaceState(null,"","#"+v);}catch(e){}
  if(v!=="map") window.scrollTo({top:0});
}
document.querySelectorAll("nav.tabs button").forEach(b=>b.addEventListener("click",()=>setView(b.dataset.view)));
setView((location.hash||"#map").slice(1));
window.addEventListener("hashchange",()=>setView((location.hash||"#map").slice(1)));
})();
