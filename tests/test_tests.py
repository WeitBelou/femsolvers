import pytest


class TestTestsSolver:
    def testSuccess(self):
        """
        Tests that tests success end
        """
        assert 0 == 0

    def testFails(self):
        """
        Tests that tests fails end
        """
        assert 0 == 1
