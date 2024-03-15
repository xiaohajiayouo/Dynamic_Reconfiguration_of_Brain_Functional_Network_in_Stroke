import os
import argparse
import numpy as np
from sklearn.metrics import silhouette_score


def main():
    parser = argparse.ArgumentParser(description='Compute the silhouette score of the latent space')
    parser.add_argument('--expdir', type=str, default='results/graph_neural_mapping', help='path to the experiment results')
    parser.add_argument('--latentdir', type=str, default='latent', help='path containing the latent_space_*.npy')
    parser.add_argument('--savedir', type=str, default='silhouette', help='path to save the silhouette value within the expdir')
    parser.add_argument('--fold_idx', nargs='+', default=['0','1','2','3','4','5','6','7','8','9'], help='fold indices')

    opt = parser.parse_args()
    os.makedirs(os.path.join(opt.expdir, opt.savedir), exist_ok=True)

    latent_space_initial = []
    latent_space = []
    labels = []

    for current_fold in opt.fold_idx:
        latent_space_initial.append(np.load(os.path.join(opt.expdir, opt.latentdir, str(current_fold), 'latent_space_initial.npy')))
        latent_space.append(np.load(os.path.join(opt.expdir, opt.latentdir, str(current_fold), 'latent_space.npy')))
        labels.append(np.load(os.path.join(opt.expdir, opt.latentdir, str(current_fold), 'labels.npy')).squeeze())

    initial_silhouette = []
    silhouette = []

    for idx, _ in enumerate(labels):
        initial_silhouette.append(silhouette_score(latent_space_initial[idx], labels[idx], sample_size=latent_space_initial[idx].shape[0]))
        silhouette.append(silhouette_score(latent_space[idx], labels[idx], sample_size=latent_space[idx].shape[0]))

    with open(os.path.join(opt.expdir, opt.savedir,'silhouette_score.csv'), 'w') as f:
        f.write("fold_idx,initial_silhouette_score,silhouette_score\n")
        for idx, (init_sil, sil) in enumerate(zip(initial_silhouette, silhouette)):
            f.write("{},{},{}\n".format(idx, init_sil, sil))
        f.write("mean,{},{}\n".format(np.mean(initial_silhouette), np.mean(silhouette)))


if __name__=='__main__':
    main()
