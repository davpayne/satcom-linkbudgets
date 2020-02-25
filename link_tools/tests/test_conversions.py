from unittest import TestCase

import link_tools

class TestERPConvert(TestCase):
    def eirp_to_erp_convert(self):
        s = link_tools.EIRP_to_ERP_dBW(2.15)
        self.assertTrue(isinstance(s, 0))
