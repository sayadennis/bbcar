# pylint: disable=missing-module-docstring
# pylint: disable=too-many-arguments

import matplotlib.pyplot as plt
from suphNMF import suphNMF


def plot_learning(
    X1,
    X2,
    y,
    n_iter: int,
    lr: float,
    weight_decay: float,
    clf_weight: float,
    ortho_weight: float,
    plotdir: str,
):
    """
    Plot learning curves for the respective loss elements.

    Parameters:
        weight_decay
        clf_weight
        ortho_weight

    Returns:
        None
    """
    test = suphNMF(X1, X2, y)
    test.fit(
        X1,
        X2,
        y,
        n_iter=n_iter,
        lr=lr,
        weight_decay=weight_decay,
        clf_weight=clf_weight,
        ortho_weight=ortho_weight,
    )

    fig, axs = plt.subplots(2, 3, figsize=(12, 8))

    axs[0, 0].plot(test.loss_record)
    axs[0, 0].set_title("Overall Loss")

    axs[0, 1].plot(test.reconloss1_record)
    axs[0, 1].set_title("Reconstruction Loss (Mutation)")

    axs[0, 2].plot(test.reconloss2_record)
    axs[0, 2].set_title("Reconstruction Loss (CNV)")

    axs[1, 0].plot(test.ortholoss_record)
    axs[1, 0].set_title("Orthogonality loss")

    axs[1, 1].plot(test.clf_roc_record)
    axs[1, 1].set_ylabel("Classification performance")
    axs[1, 1].set_title("Classification ROC-AUC")

    axs[1, 2].plot(test.clf_loss_record)
    axs[1, 2].set_title("Classification Loss")

    for i in range(2):
        for j in range(3):
            axs[i, j].set_xlabel("Epochs")
            if not (i == 1) & (j == 1):
                axs[i, j].set_ylabel("Loss")

    fig.suptitle(
        f"Learning curves - weight decay {weight_decay}; "
        f"classification weight {clf_weight}; "
        f"orthogonality weight {ortho_weight}"
    )
    plt.tight_layout()
    fig.savefig(
        f"{plotdir}/test_suphNMF_learning_curves_decay{weight_decay}"
        f"_clf{clf_weight}_ortho{ortho_weight}.png"
    )
    plt.close()
