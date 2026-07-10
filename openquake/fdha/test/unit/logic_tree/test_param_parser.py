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


def test_parse_oq_engine_style_element_text():
    # Plain element text as in oq-engine GMPE logic trees: block indented to
    # the XML nesting depth, no CDATA. Irregular (even increasing) per-line
    # indentation must not turn a "key = value" line into a continuation of
    # the previous value.
    text = (
        "\n            [Youngs2003SecondarySR]"
        "\n            version = 3"
        "\n                style = all"
        "\n          "
    )
    cls, params = parse_uncertainty_model(text)
    assert cls == "Youngs2003SecondarySR"
    assert params == {"version": 3, "style": "all"}


def test_parse_plain_class_name_padded():
    cls, params = parse_uncertainty_model("\n        Chiou2025PrimaryFD\n    ")
    assert cls == "Chiou2025PrimaryFD"
    assert params == {}

