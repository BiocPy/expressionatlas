import pytest
import rds2py
from expressionatlas.download import _download_and_load_rds

def test_download_and_load_r_files(tmp_path):
    # Test rds loading
    rds_path = tmp_path / "test.rds"
    rds2py.write_rds({"value": [1, 2, 3]}, str(rds_path))
    
    # Test loading via local file URL
    file_url = rds_path.as_uri()
    res = _download_and_load_rds(file_url, "dummy")
    assert res.names is not None and "value" in list(res.names)
    assert [list(x)[0] for x in res["value"]] == [1, 2, 3]

    # Test rda loading
    rdata_path = tmp_path / "test.Rdata"
    rds2py.write_rda({"value": [4, 5, 6]}, str(rdata_path))
    
    file_url_rdata = rdata_path.as_uri()
    res_rdata = _download_and_load_rds(file_url_rdata, "dummy")
    assert res_rdata.names is not None and "value" in list(res_rdata.names)
    assert [list(x)[0] for x in res_rdata["value"]] == [4, 5, 6]
