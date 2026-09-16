"""Toy evidence-bundle tests. No production archive or numerical replay."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import zipfile
import imt_reproduce as imt


class IMTBundleTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='imt-bundle-test-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        source = self.root/'source.txt'; source.write_text('test input',encoding='utf-8')
        record = {'source_sha256': {'source.txt':imt.sha(source)}}
        receipt = self.root/'receipt.json'; receipt.write_text(json.dumps(record),encoding='utf-8')
        self.roots = {'receipt.json':imt.sha(receipt)}

    def test_inventory_and_archive(self):
        report = imt.inventory(self.root,self.roots)
        self.assertTrue(report['inventory_complete'])
        self.assertEqual((report['entries'],report['missing'],report['mismatched']),(2,0,0))
        result = imt.pack(self.root/'evidence.zip',self.root,self.roots)
        self.assertEqual(result['entries'],2)
        imt.verify_archive(self.root/'evidence.zip',report)

    def test_missing_root_is_explicit(self):
        report = imt.inventory(self.root,{'missing.json':'0'*64})
        self.assertFalse(report['inventory_complete'])
        self.assertEqual(report['unreadable_roots'],['missing.json'])
        with self.assertRaises(ValueError):
            imt.pack(self.root/'bad.zip',self.root,{'missing.json':'0'*64})

    def test_changed_dependency_rejected(self):
        (self.root/'source.txt').write_text('changed')
        self.assertEqual(imt.inventory(self.root,self.roots)['mismatched'],1)
        with self.assertRaises(ValueError):
            imt.pack(self.root/'bad.zip',self.root,self.roots)

    def test_never_overwrite(self):
        output = self.root/'existing.zip'; output.write_bytes(b'preserve')
        with self.assertRaises(FileExistsError):
            imt.pack(output,self.root,self.roots)
        self.assertEqual(output.read_bytes(),b'preserve')

    def test_paths_and_conflicting_pins(self):
        with self.assertRaises(ValueError): imt.checked_path(self.root,'../escape')
        with self.assertRaises(ValueError): imt.checked_path(self.root,'bad\\path')
        with self.assertRaises(ValueError): imt.add_pin({'a':'0'*64},'a','1'*64)

    def test_archive_tamper_detected(self):
        output = self.root/'tampered.zip'
        report = imt.inventory(self.root,self.roots)
        with zipfile.ZipFile(output,'w') as archive:
            archive.writestr('source.txt',b'wrong')
            archive.writestr('receipt.json',(self.root/'receipt.json').read_bytes())
            archive.writestr('imt-evidence-manifest.json',json.dumps(report))
        with self.assertRaises(ValueError): imt.verify_archive(output,report)


if __name__ == '__main__':
    unittest.main()
