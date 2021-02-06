import os
from ligbinder.settings import Settings
from .conftest import FPATH

def test_settings_defaults():
    s = Settings()
    assert s["tree"]["max_depth"] == 20
    assert len(s) == 5
    assert "md" in s.data

def test_settings_override():
    file1 = os.path.join(FPATH, "settings", "settings1.yml") 
    s = Settings(file1)
    assert s["my_setting"] == "test_setting"
    assert s["md"]["crd_file"] == "test_crd_file.crd"
    assert len(s["md"]) > 1
    assert s["md"]["top_file"] == "top.prmtop"