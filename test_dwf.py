import digitalwellenfilter as dwf

def test_Element1_init():
    e = dwf.Element(1, 1)
    assert e.a == 1
    assert e.delay[0] == 0
    assert e.delay[1] == 0
    assert e.ein == 0
    assert e.aus == 0
    assert e.typ == 1


def test_Element1_advance():
    a = 0.333
    e = dwf.Element(a, 1)
    e.ein = 1
    e.advance()
    assert e.aus == a
    assert e.delay[0] == 1.333
    assert e.delay[1] == 0

    e.ein = 0
    e.advance()
    assert e.aus == 0
    assert e.delay[0] == 0
    assert e.delay[1] == 1.333

    e.advance()
    assert e.aus == 1.333 * a - 1.333
    assert e.delay[0] == 1.333 + a * 1.333 - 1.333
    assert e.delay[1] == 0


def test_Element2_advance():
    a = 0.333
    e = dwf.Element(a, 2)
    e.ein = 1
    e.advance()
    assert e.aus == -1 + a * 1
    assert e.delay[0] == a * 1
    assert e.delay[1] == 0

    e.ein = 0
    e.advance()
    assert e.aus == 0
    assert e.delay[0] == 0
    assert e.delay[1] == a * 1

    e.advance()
    assert e.aus == -a + (-1 + a) * a
    assert e.delay[0] == -a + a*a
    assert e.delay[1] == 0
