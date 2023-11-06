"""Utility functions for sanity checks."""


def run_incompatible_value_checks(
    pseudo_bulk, loss_computation, use_batch_norm, mixup_penalty, gene_likelihood
):
    """Check the values of the categorical variables to run MixUpVI are compatible.
    The first 4 checks will only be relevant when pseudobulk will not be computed both
    in encoder and decoder (right now, computed in both). Until then, use_batch_norm
    should be None.
    """
    if (
        pseudo_bulk == "pre_encoded"
        and loss_computation == "latent_space"
        and use_batch_norm in ["encoder", "both"]
    ):
        raise ValueError(
            "MixUpVI cannot use batch normalization there, as the batch size of pseudobulk is 1."
        )
    elif (
        pseudo_bulk == "pre_encoded"
        and loss_computation == "reconstructed_space"
        and use_batch_norm != "none"
    ):
        raise ValueError(
            "MixUpVI cannot use batch normalization there, as the batch size of pseudobulk is 1."
        )
    elif pseudo_bulk == "post_inference" and loss_computation == "latent_space":
        raise ValueError(
            "Pseudo bulk needs to be pre-encoded to compute the MixUp loss in the latent space."
        )
    elif (
        pseudo_bulk == "post_inference"
        and loss_computation == "reconstructed_space"
        and use_batch_norm in ["decoder", "both"]
    ):
        raise ValueError(
            "MixUpVI cannot use batch normalization there, as the batch size of pseudobulk is 1."
        )
    if (
        mixup_penalty == "kl"
        and loss_computation != "latent_space"
        and gene_likelihood == "zinb"
    ):
        raise NotImplementedError(
            "The KL divergence between ZINB distributions for the MixUp loss is not "
            "implemented."
        )


def run_categorical_value_checks(
    cell_group,
    cat_cov,
    cont_cov,
    encode_covariates,
    encode_cont_covariates,
    use_batch_norm,
    signature_type,
    loss_computation,
    pseudo_bulk,
    mixup_penalty,
    dispersion,
    gene_likelihood,
):
    """Check the values and types of the categorical variables to run MixUpVI."""
    assert isinstance(cell_group, str), "CELL_GROUP should be of type string"
    assert isinstance(cat_cov, list), "CAT_COV should be of type list"
    assert (
        isinstance(cont_cov, list) or cont_cov == None
    ), "CONT_COV should be None or type list"
    assert isinstance(
        encode_covariates, bool
    ), "ENCODE_COVARIATES should be of type bool"
    assert isinstance(
        encode_cont_covariates, bool
    ), "ENCODE_CONT_COVARIATES should be of type bool"
    assert isinstance(use_batch_norm, str), "BATCH_NORM should be of type string"
    assert isinstance(signature_type, str), "SIGNATURE_TYPE should be of type string"
    assert isinstance(
        loss_computation, str
    ), "LOSS_COMPUTATION should be of type string"
    assert isinstance(pseudo_bulk, str), "PSEUDO_BULK should be of type string"
    assert isinstance(mixup_penalty, str), "MIXUP_PENALTY should be of type string"
    assert isinstance(dispersion, str), "DISPERSION should be of type string"
    assert isinstance(gene_likelihood, str), "GENE_LIKELIHOOD should be of type string"
    if cell_group not in [
        "primary_groups",
        "precise_groups",
        "updated_granular_groups",
    ]:
        raise NotImplementedError(
            "For now, the following cell category granularities are implemented: "
            "['primary_groups', 'precise_groups', 'updated_granular_groups']"
        )
    if len(cat_cov) > 1:
        raise NotImplementedError(
            "For now, MixUpVI works with only the cell type as categorical covariate."
        )
    if encode_covariates:
        raise NotImplementedError(
            "For now, MixUpVI only uses cell types as categorical covariates without encoding them."
        )
    if use_batch_norm not in ["encoder", "decoder", "none", "both"]:
        raise ValueError(
            "Batch normalization can only be part of ['encoder', 'decoder', 'none', 'both']."
        )
    if signature_type not in ["pre_encoded", "post_inference"]:
        raise ValueError(
            "Signature type can only be part of ['pre_encoded', 'post_inference']."
        )
    if loss_computation not in ["latent_space", "reconstructed_space"]:
        raise ValueError(
            "Loss computation can only be part of ['latent_space', 'reconstructed_space']."
        )
    if pseudo_bulk not in ["pre_encoded", "post_inference"]:
        raise ValueError(
            "Pseudo bulk computation can only be part of ['pre_encoded', 'post_inference']."
        )
    if mixup_penalty not in ["l2", "kl"]:
        raise ValueError("Mixup penalty can only be part of ['l2', 'kl'].")
    if dispersion not in ["gene", "gene_cell"]:
        raise ValueError(
            "The dispersion parameter can only be part of ['gene', 'gene_cell'], "
            "not gene-label nor gene-batch because categorical covariates don't make "
            "sense for pseudobulk."
        )
    if gene_likelihood not in ["zinb", "nb", "poisson"]:
        raise ValueError(
            "The dispersion parameter can only be part of ['zinb', 'nb', 'poisson']."
        )