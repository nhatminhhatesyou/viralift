from app.src.alias.alias_classifier import classify_alias_candidate


def _classify(raw_value: str, canonical_name: str):
    return classify_alias_candidate(
        raw_value=raw_value,
        field="product",
        canonical_name=canonical_name,
        evidence={"iou": 0.96, "strand_match": True},
    )


def _classify_field(raw_value: str, canonical_name: str, field: str):
    return classify_alias_candidate(
        raw_value=raw_value,
        field=field,
        canonical_name=canonical_name,
        evidence={"iou": 0.96, "strand_match": True},
    )


def test_orf1a_protein_matches_orf1a():
    result = _classify("ORF1a protein", "ORF1a")

    assert result["suggested_action"] == "save_alias"
    assert result["confidence"] == "high"
    assert "contains canonical name" in result["reason"]


def test_orf1a_protein_does_not_match_orf1ab_by_short_symbol():
    result = _classify("ORF1a protein", "ORF1ab")

    assert "short gene symbol consistent with canonical" not in result["reason"]
    assert result["confidence"] != "high"


def test_envelope_protein_matches_e_canonical():
    result = _classify("envelope protein", "E")

    assert result["suggested_action"] == "save_alias"
    assert result["confidence"] == "high"
    assert "descriptive synonym matches canonical" in result["reason"]


def test_membrane_protein_matches_m_canonical():
    result = _classify("membrane protein", "M")

    assert result["suggested_action"] == "save_alias"
    assert result["confidence"] == "high"
    assert "descriptive synonym matches canonical" in result["reason"]


def test_membrane_protein_does_not_match_e_canonical():
    result = _classify("membrane protein", "E")

    assert "descriptive synonym matches canonical" not in result["reason"]
    assert result["suggested_action"] != "save_alias"


def test_contextual_descriptions_are_review_terms_not_hard_blacklist():
    """Names whose meaning depends on context stay reviewable, not blacklisted.

    'glycoprotein' is ambiguous rather than unusable: with an anchoring number it
    becomes a real alias, so it must reach review instead of being ruled out.
    """
    for raw_value in [
        "glycoprotein",
        "major glycoprotein",
        "minor glycoprotein",
    ]:
        result = _classify(raw_value, "ORF1a")
        assert "generic name" not in result["reason"]
        assert "descriptive biological term" in result["reason"]


def test_class_words_and_bare_activities_are_blacklisted():
    """Class words and bare activities are ruled out, not left to per-run judgement.

    These previously sat in the review bucket. That made the outcome depend on how
    the LLM felt about them: across four PRRSV reconstruction runs every false save
    came from exactly this pool, 2 to 5 accepted each time. They cannot name one
    gene in any virus — several genes share the class, and an activity with no gene
    designation names a domain — so the decision belongs here where it is stable.
    """
    for raw_value in [
        "structural protein",
        "nonstructural protein",
        "non-structural protein",
        "structural membrane protein",
        "polyprotein",
        "proteinase",
    ]:
        result = _classify(raw_value, "ORF1a")
        assert "generic name" in result["reason"], raw_value
        assert result["suggested_action"] == "ignore", raw_value


def test_blacklist_does_not_catch_labels_carrying_a_gene_designation():
    """The digit guard keeps numbered forms out of the blacklist."""
    for raw_value, canonical in [
        ("polyprotein 1a", "ORF1a"),
        ("non-structural protein ORF1a", "ORF1a"),
        ("structural protein GP5", "ORF5"),
    ]:
        result = _classify(raw_value, canonical)
        assert "generic name" not in result["reason"], raw_value
        assert result["suggested_action"] != "ignore", raw_value


def test_specific_polyprotein_description_matches_orf_canonical():
    cases = [
        ("polyprotein 1a", "ORF1a"),
        ("replicase polyprotein 1b", "ORF1b"),
        ("polyprotein 1ab", "ORF1ab"),
    ]

    for raw_value, canonical_name in cases:
        result = _classify(raw_value, canonical_name)
        assert result["suggested_action"] == "save_alias"
        assert result["confidence"] == "high"
        assert "specific ORF polyprotein description matches canonical" in result["reason"]


def test_specific_polyprotein_description_does_not_match_wrong_orf_canonical():
    result = _classify("polyprotein 1a", "ORF1ab")

    assert "specific ORF polyprotein description matches canonical" not in result["reason"]
    assert result["suggested_action"] != "save_alias"


def test_hard_noise_names_remain_generic():
    for raw_value in ["unknown", "hypothetical protein", "protein"]:
        result = _classify(raw_value, "ORF3")
        assert result["suggested_action"] == "ignore"
        assert result["confidence"] == "low"
        assert "generic name" in result["reason"]


def test_gene_field_without_name_specific_evidence_is_not_high_confidence():
    for raw_value in ["HNZK1", "mp"]:
        result = _classify_field(raw_value, "ORF3", "gene")
        assert result["suggested_action"] == "manual_review"
        assert result["confidence"] == "medium"
        assert "no name-specific alias evidence" in result["reason"]


def test_short_exact_gene_symbol_can_still_be_high_confidence():
    result = _classify_field("S", "S", "gene")

    assert result["suggested_action"] == "save_alias"
    assert result["confidence"] == "high"
    assert "exact canonical text" in result["reason"]
