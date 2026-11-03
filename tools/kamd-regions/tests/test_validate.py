import numpy as np

from conftest import synthetic_membrane
from kamd_regions.analysis import Analysis
from kamd_regions.project import Project, Region
from kamd_regions.system import load_system
from kamd_regions.validate import validate


def codes(v, key="errors"):
    return {x["code"] for x in v[key]}


def test_amber_import_is_valid(amber_session):
    v = amber_session.validate()
    assert not v["errors"]
    assert "small_group" in codes(v, "warnings")


def test_error_cases(amber_session):
    an = amber_session.an
    n = an.s.n_residues
    lab = np.zeros(n, dtype=int)
    lab[1:5] = 1
    p = Project.from_labels(lab)
    assert "too_few_groups" in codes(validate(an, p))
    p.regions.append(Region(9, "empty", "#000000"))
    assert "empty_group" in codes(validate(an, p))
    # > 32 groups is an Amber-only error
    lab = np.zeros(n, dtype=int)
    lab[:39] = np.arange(1, 40)
    p = Project.from_labels(lab)
    v = validate(an, p)
    assert "too_many_groups" in codes(v)
    assert v["ok"] == {"amber": False, "gromacs": True, "openmm": True}
    # transmembrane group without membrane
    lab = np.zeros(n, dtype=int)
    lab[1:10] = 1
    lab[20:30] = 2
    p = Project.from_labels(lab, {1: 2, 2: 1})
    assert {"no_membrane_bath", "no_membrane_axis"} <= codes(validate(an, p))
    # water in a group
    lab[100] = 2
    assert "bath_in_group" in codes(validate(an, Project.from_labels(lab)))


def test_links(amber_session):
    p = amber_session.project
    q = Project.from_dict(p.to_dict())
    for r in q.regions[1:]:
        r.refer_link = q.regions[0].id
    assert "all_linked" in codes(validate(amber_session.an, q), "warnings")
    q.regions[0].refer_link = q.regions[0].id
    assert "self_link" in codes(validate(amber_session.an, q))


def test_membrane_warnings_and_contiguity(tmp_path):
    s = load_system([str(synthetic_membrane(tmp_path, "lipid21", split_water=True))])
    an = Analysis(s)
    lab = np.zeros(s.n_residues, dtype=int)
    lab[0] = 1
    lab[1:3] = 2
    p = Project.from_labels(lab, {1: 1, 2: 2})
    v = validate(an, p)
    assert "solvent_not_contiguous" in codes(v)
    assert v["ok"]["amber"] is False and v["ok"]["openmm"] is True
    assert "no_topology" in codes(v, "warnings")
