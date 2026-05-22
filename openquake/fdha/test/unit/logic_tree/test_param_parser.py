from openquake.fdha.logic_tree.param_parser import parse_uncertainty_model


def test_parse_plain_class_name():
    cls, params = parse_uncertainty_model("Chiou2025PrimaryFD")
    assert cls == "Chiou2025PrimaryFD"
    assert params == {}


def test_parse_ini_block_literals():
    text = """
    [Pizza2023PrimarySR]
    style = "all"
    MSR = 1
    """
    cls, params = parse_uncertainty_model(text)
    assert cls == "Pizza2023PrimarySR"
    assert params["style"] == "all"
    assert params["MSR"] == 1

