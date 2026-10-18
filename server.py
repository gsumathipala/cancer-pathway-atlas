#!/usr/bin/env python3
"""Cancer Pathway Atlas: a dependency-free web service.

Serves the interactive atlas from ./docs and a read-only JSON API built from docs/data.json.

    python3 server.py                  # http://127.0.0.1:8000
    python3 server.py --port 9000 --open
    PORT=8080 python3 server.py --host 0.0.0.0
"""
import argparse
import json
import os
import re
import sys
import webbrowser
from functools import partial
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

ROOT = Path(__file__).resolve().parent / "docs"
EDGE_RE = re.compile(r"^(\w+)([>|\-])(\w+)$")
EDGE_KIND = {">": "activates", "|": "inhibits", "-": "complex"}
STATUS = {"A": "approved", "T": "trials", "X": "stopped"}
CSP = ("default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline' https://fonts.googleapis.com; "
       "font-src https://fonts.gstatic.com; img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'")


class Atlas:
    """In-memory index over docs/data.json."""

    def __init__(self, path):
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        self.meta = raw["meta"]
        self.refs = [{"n": i + 1, "citation": c, "doi": d, "url": "https://doi.org/" + d} for i, (c, d) in enumerate(raw["refs"])]
        self.order = raw["order"]["tcga"] + raw["order"]["other"]
        self.tcga = set(raw["order"]["tcga"])
        keys = ["id", "label", "pathway", "role", "frequency", "x", "y", "alterations", "note", "cancers", "drugs"]
        self.genes = {g["id"]: g for g in (dict(zip(keys, row)) for row in raw["nodes"])}
        self.edges = []
        for token in raw["edges"].split():
            m = EDGE_RE.match(token)
            if m:
                self.edges.append({"source": m.group(1), "target": m.group(3), "type": EDGE_KIND[m.group(2)]})
        self.pathways = {}
        for pid, p in raw["pathways"].items():
            self.pathways[pid] = {
                "id": pid, "name": p["name"], "short": p["short"], "canonical": pid in self.tcga,
                "frequency": p["fr"], "summary": p["mech"], "lesions": p["lesion"], "where": p["where"],
                "vulnerabilities": p["vuln"], "references": p["refs"],
                "steps": [{"title": s[0], "text": s[1], "refs": s[2]} for s in p.get("steps", [])],
                "evidence": [{"year": e[0], "finding": e[1], "refs": e[2]} for e in p.get("evidence", [])],
                "therapies": [self._therapy(pid, t) for t in p.get("therapies", []) if t[1] != "—"],
            }
        self.by_label = {}
        for g in self.genes.values():
            self.by_label[g["id"].upper()] = g["id"]
            self.by_label[g["label"].upper()] = g["id"]
            for part in re.split(r"[\s/()–-]+", g["label"].upper()):
                if part:
                    self.by_label.setdefault(part, g["id"])

    @staticmethod
    def _therapy(pid, t):
        return {"pathway": pid, "agent": t[0], "target": t[1], "setting": t[2], "status": STATUS[t[3]],
                "fdaYear": int(t[4]) if t[4] else None, "refs": t[5]}

    def pathway_summary(self, pid):
        p = self.pathways[pid]
        return {k: p[k] for k in ("id", "name", "short", "canonical", "frequency")} | {
            "genes": sum(1 for g in self.genes.values() if g["pathway"] == pid),
            "approvedTherapies": sum(1 for t in p["therapies"] if t["status"] == "approved")}

    def resolve_gene(self, key):
        key = key.upper()
        return self.by_label.get(key)

    def gene_detail(self, gid):
        g = dict(self.genes[gid])
        g.pop("x", None)
        g.pop("y", None)
        g["partners"] = [
            {"gene": e["target"] if e["source"] == gid else e["source"], "type": e["type"],
             "direction": "out" if e["source"] == gid else "in"}
            for e in self.edges if gid in (e["source"], e["target"])]
        return g

    def therapies(self, pathway=None, status=None):
        rows = [t for pid in self.order for t in self.pathways[pid]["therapies"]]
        if pathway:
            rows = [t for t in rows if t["pathway"] == pathway]
        if status:
            rows = [t for t in rows if t["status"] == status]
        return rows

    def search(self, q):
        q = q.strip().lower()
        if not q:
            return {"genes": [], "pathways": [], "therapies": []}
        genes = [{"id": g["id"], "label": g["label"], "pathway": g["pathway"]} for g in self.genes.values()
                 if q in g["label"].lower() or q in g["id"].lower() or q in g["cancers"].lower()]
        paths = [self.pathway_summary(pid) for pid, p in self.pathways.items()
                 if q in p["name"].lower() or q in p["short"].lower()]
        ther = [t for t in self.therapies() if q in (t["agent"] + " " + t["target"] + " " + t["setting"]).lower()]
        return {"genes": genes, "pathways": paths, "therapies": ther}


class Handler(SimpleHTTPRequestHandler):
    atlas = None
    server_version = "CancerPathwayAtlas/1.1"

    def end_headers(self):
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy", CSP)
        super().end_headers()

    def log_message(self, fmt, *args):
        sys.stderr.write("%s %s\n" % (self.address_string(), fmt % args))

    def do_GET(self):
        url = urlparse(self.path)
        if url.path == "/api" or url.path.startswith("/api/"):
            return self.api(url)
        return super().do_GET()

    def do_HEAD(self):
        if urlparse(self.path).path.startswith("/api"):
            return self.send_json({"error": "use GET"}, HTTPStatus.METHOD_NOT_ALLOWED)
        return super().do_HEAD()

    def send_json(self, obj, status=HTTPStatus.OK):
        body = json.dumps(obj, ensure_ascii=False, indent=1).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(body)

    def api(self, url):
        a = self.atlas
        parts = [unquote(p) for p in url.path.split("/")[2:] if p]
        qs = {k: v[0] for k, v in parse_qs(url.query).items()}
        notfound = lambda what: self.send_json({"error": "%s not found" % what}, HTTPStatus.NOT_FOUND)
        if not parts:
            return self.send_json({"endpoints": [
                "/api/health", "/api/pathways", "/api/pathways/{id}", "/api/genes?pathway={id}",
                "/api/genes/{symbol}", "/api/edges", "/api/therapies?pathway={id}&status=approved|trials|stopped",
                "/api/search?q={text}", "/api/references"]})
        head, rest = parts[0], parts[1:]
        if head == "health" and not rest:
            return self.send_json({"status": "ok", "version": a.meta["version"], "pathways": len(a.pathways),
                                   "genes": len(a.genes), "edges": len(a.edges), "references": len(a.refs)})
        if head == "pathways":
            if not rest:
                return self.send_json([a.pathway_summary(pid) for pid in a.order])
            pid = rest[0].lower()
            if pid not in a.pathways:
                return notfound("pathway")
            genes = [a.gene_detail(g["id"]) for g in a.genes.values() if g["pathway"] == pid]
            return self.send_json(a.pathways[pid] | {"genes": genes})
        if head == "genes":
            if not rest:
                pw = qs.get("pathway")
                rows = [g for g in a.genes.values() if not pw or g["pathway"] == pw]
                return self.send_json([{k: g[k] for k in ("id", "label", "pathway", "role", "frequency")} for g in rows])
            gid = a.resolve_gene(rest[0])
            return self.send_json(a.gene_detail(gid)) if gid else notfound("gene")
        if head == "edges" and not rest:
            return self.send_json(a.edges)
        if head == "therapies" and not rest:
            status = qs.get("status")
            if status and status not in STATUS.values():
                return self.send_json({"error": "status must be approved, trials or stopped"}, HTTPStatus.BAD_REQUEST)
            return self.send_json(a.therapies(qs.get("pathway"), status))
        if head == "search" and not rest:
            return self.send_json(a.search(qs.get("q", "")))
        if head == "references" and not rest:
            return self.send_json(a.refs)
        return notfound("endpoint")


def make_server(host, port, root=ROOT):
    Handler.atlas = Atlas(Path(root) / "data.json")
    return ThreadingHTTPServer((host, port), partial(Handler, directory=str(root)))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Serve the Cancer Pathway Atlas.")
    ap.add_argument("--host", default=os.environ.get("HOST", "127.0.0.1"))
    ap.add_argument("--port", type=int, default=int(os.environ.get("PORT", "8000")))
    ap.add_argument("--open", action="store_true", help="open the atlas in a browser")
    args = ap.parse_args(argv)
    httpd = make_server(args.host, args.port)
    url = "http://%s:%d/" % ("127.0.0.1" if args.host in ("0.0.0.0", "") else args.host, args.port)
    print("Cancer Pathway Atlas running at %s  (API: %sapi)  Ctrl+C to stop" % (url, url))
    if args.open:
        webbrowser.open(url)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()


if __name__ == "__main__":
    main()
