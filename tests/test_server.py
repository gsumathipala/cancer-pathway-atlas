import json
import sys
import threading
import unittest
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import server  # noqa: E402


class ApiTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.httpd = server.make_server("127.0.0.1", 0)
        cls.base = "http://127.0.0.1:%d" % cls.httpd.server_address[1]
        threading.Thread(target=cls.httpd.serve_forever, daemon=True).start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()
        cls.httpd.server_close()

    def get(self, path):
        with urlopen(self.base + path) as r:
            return r.status, r.headers, r.read()

    def json(self, path):
        return json.loads(self.get(path)[2])

    def test_health(self):
        h = self.json("/api/health")
        self.assertEqual(h["status"], "ok")
        self.assertEqual(h["pathways"], 17)
        self.assertGreater(h["genes"], 90)

    def test_pathway_detail_has_tutorial_and_therapies(self):
        for p in self.json("/api/pathways"):
            d = self.json("/api/pathways/" + p["id"])
            self.assertTrue(d["steps"], p["id"])
            self.assertTrue(d["evidence"], p["id"])
            self.assertTrue(d["genes"], p["id"])

    def test_all_reference_numbers_resolve(self):
        n = len(self.json("/api/references"))
        for p in self.json("/api/pathways"):
            d = self.json("/api/pathways/" + p["id"])
            nums = d["references"] + [r for s in d["steps"] for r in s["refs"]] + [r for e in d["evidence"] for r in e["refs"]] \
                + [r for t in d["therapies"] for r in t["refs"]]
            self.assertTrue(all(1 <= r <= n for r in nums), p["id"])

    def test_gene_lookup_by_symbol_and_alias(self):
        self.assertEqual(self.json("/api/genes/tp53")["pathway"], "p53")
        self.assertEqual(self.json("/api/genes/HER2")["id"], "ERBB2")
        self.assertTrue(any(x["gene"] == "MDM2" for x in self.json("/api/genes/TP53")["partners"]))

    def test_therapy_filters(self):
        appr = self.json("/api/therapies?status=approved")
        self.assertTrue(appr and all(t["status"] == "approved" for t in appr))
        self.assertTrue(all(t["pathway"] == "ddr" for t in self.json("/api/therapies?pathway=ddr")))
        with self.assertRaises(HTTPError) as cm:
            self.get("/api/therapies?status=bogus")
        self.assertEqual(cm.exception.code, 400)

    def test_unknown_routes_404(self):
        for path in ("/api/genes/NOTAGENE", "/api/pathways/nope", "/api/nothing"):
            with self.assertRaises(HTTPError) as cm:
                self.get(path)
            self.assertEqual(cm.exception.code, 404)

    def test_static_site_and_headers(self):
        status, headers, body = self.get("/")
        self.assertEqual(status, 200)
        self.assertIn(b"Cancer Pathway Atlas", body)
        self.assertIn("default-src 'self'", headers["Content-Security-Policy"])
        self.assertEqual(self.get("/vendor/d3.min.js")[0], 200)

    def test_no_path_traversal(self):
        with self.assertRaises(HTTPError) as cm:
            self.get("/../server.py")
        self.assertEqual(cm.exception.code, 404)


if __name__ == "__main__":
    unittest.main()
