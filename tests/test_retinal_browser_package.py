import json
from pathlib import Path
import sys
import tempfile
import tarfile
import unittest
import pyBigWig
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from assemble_retinal_browser_package import verify_group, make_archive
from prepare_browser_tracks import digest


class RetinalBrowserPackage(unittest.TestCase):
    def test_archive_preserves_content_without_machine_owner_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);folder=root/'review';folder.mkdir()
            (folder/'signal.bw').write_bytes(b'example track content')
            archive=root/'review.tar'
            metadata=make_archive(folder,archive)
            self.assertEqual(metadata['sha256'],digest(archive))
            with tarfile.open(archive) as tar:
                entry=tar.getmember('review/signal.bw')
                self.assertEqual((entry.uid,entry.gid,entry.uname,entry.gname),(0,0,'',''))
                self.assertEqual(tar.extractfile(entry).read(),b'example track content')
            with self.assertRaises(FileExistsError):make_archive(folder,archive)

    def fixture(self,root):
        directory=root/'MG';directory.mkdir()
        donors={d:{'fragments':1,'cells':1,'sha256':'input-'+d} for d in ['LGS1','LGS2','LGS3','LVG1']}
        config={'groups':{'MG':{'donors':donors,'keep_donor_tracks':True}}}
        records=[]
        for donor in ['pooled']+list(donors):
            path=directory/('MG_'+donor+'.IPM_per_base.100bp.bw')
            with pyBigWig.open(str(path),'w') as bw:
                bw.addHeader([('chr1',100)])
                bw.addEntries(['chr1'],[0],ends=[100],values=[10000.])
            records.append({'file':path.name,'sha256':digest(path),'insertions':8 if donor=='pooled' else 2,'bytes':path.stat().st_size})
        report={'status':'OBSERVED_CELLTYPE_BROWSER_TRACKS_READY','group':'MG','assembly':'hg38','bin_bp':100,'config_sha256':'config-hash','donors':[dict(donor=d,**v) for d,v in donors.items()],
            'fragments':4,'insertions':8,'mg_reference_check':'PASS_EXACT_ALL_GENOME_INTEGER_BIN_MATCH','tracks':records}
        self.save(directory,report)
        return directory,config,report

    def save(self,directory,report):
        (directory/'browser_track_summary.json').write_text(json.dumps(report))
        (directory/'output_checksums.sha256').write_text(''.join(digest(p)+'  ./'+p.name+'\n' for p in sorted(directory.iterdir()) if p.name!='output_checksums.sha256'))

    def test_valid_inventory_and_mg_reference_required(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);group,config,report=self.fixture(root)
            _,tracks=verify_group(root,'MG',config,'config-hash',{'chr1':100})
            self.assertEqual(len(tracks),5)
            report['mg_reference_check']=None;self.save(group,report)
            with self.assertRaises(ValueError):verify_group(root,'MG',config,'config-hash',{'chr1':100})

    def test_donor_provenance_and_file_corruption_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);group,config,report=self.fixture(root)
            report['donors'][0]['fragments']=2;self.save(group,report)
            with self.assertRaises(ValueError):verify_group(root,'MG',config,'config-hash',{'chr1':100})
            report['donors'][0]['fragments']=1;self.save(group,report)
            path=group/report['tracks'][0]['file']
            with path.open('ab') as stream:stream.write(b'corrupted')
            with self.assertRaises(ValueError):verify_group(root,'MG',config,'config-hash',{'chr1':100})

if __name__=='__main__':unittest.main()
