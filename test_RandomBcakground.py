from collections import OrderedDict
import pytest
import RandomBcakground as rb


def test_make_region_dictionary():
    expected = OrderedDict([('cds', 'CDS'),
             ('three_prime_utrs', "3' UTR"),
             ('five_prime_utrs', "5' UTR"),
             ('proxintron500', 'Proximal\nIntron'),
             ('distintron500', 'Distal\nIntron')])
    assert expected == rb.make_region_dictionary()
    
def test_check_RBP_is_a_string_1():
    rb.check_RBP_is_a_string("test")

@pytest.mark.xfail(raises=ValueError)    
def test_check_RBP_is_a_string_2():
    rb.check_RBP_is_a_string(155)
    

if __name__ == "__main__":
    pytest.main()
    
    