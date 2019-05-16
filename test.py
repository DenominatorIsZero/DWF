import digitalwellenfilter as dwf


def test_Element1_init():
    e1 = dwf.Element1(1)
    assert e1.a == 1
    assert e1.delay[0, 1] == 0, 0
    assert e1.ein == 0
    assert e1.aus == 0
