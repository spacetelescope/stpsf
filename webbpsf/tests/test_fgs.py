import logging

from .test_webbpsf import do_test_set_position_from_siaf, do_test_source_offset, generic_output_test

_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())


# ------------------    FGS Tests    ----------------------------
def test_fgs():
    return generic_output_test('FGS')


def test_fgs_source_offset_00():
    return do_test_source_offset('FGS', theta=0.0, monochromatic=2.5e-6)


def test_fgs_source_offset_45():
    return do_test_source_offset('FGS', theta=45.0, monochromatic=2.5e-6)


def test_fgs_set_siaf():
    return do_test_set_position_from_siaf('FGS',
                                          ['FGS1_FP1MIMF', 'FGS2_SUB128CNTR', 'FGS1_SUB128LL', 'FGS2_SUB32DIAG'])
