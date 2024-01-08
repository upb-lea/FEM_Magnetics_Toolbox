import numpy as np
import kalepy as kale
import pyswarms as ps
import matplotlib.pyplot as plt


class dessca_model:
    def __init__(self,
                 box_constraints,
                 bandwidth=None,
                 reference_pdf=None,
                 render_online=False,
                 PSO_options=None,
                 state_names=None,
                 buffer_size=None):

        # one time initialization
        if bandwidth is None:
            # a bandwidth of 0.1 yielded empirically good results
            # the kalepy KDE package seems to rescale this bandwidth automatically to the given box constraints and dimension
            self.bandwidth = 0.1
        else:
            self.bandwidth = bandwidth  # set the bandwidth for this estimator
        self.render_online = render_online
        self.dim = len(box_constraints)
        self.lower_bound = np.array(box_constraints)[:, 0]
        self.upper_bound = np.array(box_constraints)[:, -1]

        if PSO_options is None:
            PSO_options = {'c1': 2, 'c2': 2, 'w': 0.6}
        if state_names is None:
            self.state_names = [f"$x_{i}$" for i in range(self.dim)]
        else:
            self.state_names = state_names

        # instantiate a particle swarm optimizer
        self.optimizer = ps.single.GlobalBestPSO(n_particles=self.dim * 10,
                                                 dimensions=self.dim,
                                                 options=PSO_options,
                                                 bounds=(self.lower_bound, self.upper_bound))

        if reference_pdf is None:
            # if no referece_pdf has been defined we assume uniform distribution on the axes
            _state_space_volume = np.product([_con[-1] - _con[0] for _con in box_constraints])
            def uniform_pdf(X):
                return 1 / _state_space_volume
            self.reference_pdf = uniform_pdf
        else:
            self.reference_pdf = reference_pdf

        # ring buffer for the collected states
        if buffer_size is None:
            self.buffer_size = np.inf
            self.buffer_idx = None
            self.coverage_data = np.empty((self.dim, 0))
        else:
            self.buffer_size = buffer_size
            self.buffer_idx = 0
            self.coverage_data = np.empty((self.dim, buffer_size)) * np.nan

        # cannot be determined before data is available
        self.nb_datapoints = 0
        self.coverage_pdf = None


    def update_coverage_pdf(self, data):
        if self.render_online:
            self.render_scatter(online_data=data)
        # append the newly acquired data
        if self.buffer_idx is None:
            self.coverage_data = np.append(self.coverage_data, np.array(data), axis=1)
            coverage_data = np.copy(self.coverage_data)
        else:
            self.coverage_data[:, self.buffer_idx] = np.reshape(data, (-1))
            self.buffer_idx = (self.buffer_idx + 1) % self.buffer_size
            first_nan_idx = np.where(np.isnan(self.coverage_data))
            if len(first_nan_idx[0]) > 0:
                coverage_data = self.coverage_data[:, 0:first_nan_idx[1][0]]
            else:
                coverage_data = np.copy(self.coverage_data)


        # unfortunately this is non-recursive, recursive KDE would be more time efficient
        if np.shape(coverage_data)[1] > np.max([self.dim, 2]):
            self.coverage_pdf = kale.KDE(dataset=coverage_data, bandwidth=self.bandwidth)
        else:
            # for small no. of samples the KDE package cannot perform a density estimation,
            # hence we work around this by duplicating and adding noise to the few available samples
            _tiled = np.tile(coverage_data, (1, np.max([self.dim, 2]) + 1))
            _noisy = _tiled + np.random.normal(0, 1, np.shape(_tiled))
            self.coverage_pdf = kale.KDE(dataset=_noisy, bandwidth=self.bandwidth)

    def sample_optimally(self):
        # check if this is the first sample
        if self.coverage_pdf is not None:
            _, self.suggested_sample = self.optimizer.optimize(
                lambda X: self.coverage_pdf.density(np.transpose(X), probability=True)[1]
                          - self.reference_pdf(np.transpose(X)),
                iters=self.dim * 10 + 5,
                verbose=False)

        else:
            # if this is the first sample, the optimal sample is solely based on the reference coverage
            _, self.suggested_sample = self.optimizer.optimize(lambda X: - self.reference_pdf(np.transpose(X)),
                                                               iters=self.dim * 10 + 5,
                                                               verbose=False)
        self.optimizer.reset()

        return self.suggested_sample

    def downsample(self, data, target_size):
        # this function samples down a large dataset while preserving the original distribution
        if self.render_online:
            print("The render_online feature is not yet available for this function.")

        self.coverage_data = np.copy(data)
        dataset_size = np.shape(self.coverage_data)[1]
        self.reference_pdf = kale.KDE(dataset=self.coverage_data, bandwidth=self.bandwidth)

        _, suggested_sample = self.optimizer.optimize(
            lambda X: - self.reference_pdf.density(np.transpose(X), probability=True)[1],
            iters=self.dim * 10 + 5,
            verbose=False)
        distances = np.linalg.norm(np.add(np.transpose([suggested_sample]), -self.coverage_data), axis=0)
        removal_idx = np.argmin(distances)
        self.coverage_data = np.delete(self.coverage_data, removal_idx, 1)
        dataset_size -= 1
        self.coverage_pdf = kale.KDE(dataset=self.coverage_data, bandwidth=self.bandwidth)
        self.optimizer.reset()


        while dataset_size > target_size:
            _, suggested_sample = self.optimizer.optimize(
                lambda X: self.reference_pdf.density(np.transpose(X), probability=True)[1]
                          - self.coverage_pdf.density(np.transpose(X), probability=True)[1],
                iters=self.dim * 10 + 5,
                verbose=False)
            distances = np.linalg.norm(np.add(np.transpose([suggested_sample]), -self.coverage_data), axis=0)
            removal_idx = np.argmin(distances)
            self.coverage_data = np.delete(self.coverage_data, removal_idx, 1)
            dataset_size -= 1
            self.coverage_pdf = kale.KDE(dataset=self.coverage_data, bandwidth=self.bandwidth)
            self.optimizer.reset()


    def update_and_sample(self, data=None):
        if data is not None:
            self.update_coverage_pdf(data=data)
        return self.sample_optimally()


    def plot_heatmap(self, resolution=100, **kwargs):
        if self.dim == 1:
            print("Heatmap plot is not available for dim < 2")
        else:
            kwargs.setdefault("cmap", "inferno")
            kwargs.setdefault("vmin", 0)
            kwargs.setdefault("aspect", "equal")
            kwargs.setdefault("origin", "lower")

            self.heatmap_fig, self.heatmap_axes = plt.subplots(self.dim, self.dim)
            self.heatmap_axes = np.reshape(self.heatmap_axes, (self.dim, self.dim))

            for i in range(self.dim):
                for j in range(self.dim):
                    if j < i:
                        _xj = np.linspace(self.lower_bound[j], self.upper_bound[j], resolution)
                        _xi = np.linspace(self.lower_bound[i], self.upper_bound[i], resolution)
                        grid = np.meshgrid(_xj, _xi, indexing="ij")
                        positions = np.vstack([_grid.ravel() for _grid in grid])
                        kde = kale.KDE(np.concatenate(([self.coverage_data[j]], [self.coverage_data[i]]), axis=0))
                        _, _disc_coverage = kde.density(positions, probability=True)
                        _disc_coverage = np.reshape(_disc_coverage, grid[0].shape)
                        img = self.heatmap_axes[i, j].imshow(np.transpose(_disc_coverage), **kwargs)
                        self.heatmap_axes[i, j].set_xticks([])
                        self.heatmap_axes[i, j].set_yticks([])
                        if i == self.dim - 1:
                            self.heatmap_axes[i, j].set_xlabel(self.state_names[j])
                        if j == 0:
                            self.heatmap_axes[i, j].set_ylabel(self.state_names[i])

            self.heatmap_fig.colorbar(img, ax=self.heatmap_axes.ravel().tolist())
            for _ax in self.heatmap_axes.flatten():
                if not _ax.images:
                    _ax.remove()
            plt.show()
            return self.heatmap_fig, self.heatmap_axes

    def render_scatter(self, online_data=None):

        if not hasattr(self, "scatter_fig"):
            self.scatter_fig, self.scatter_axes = plt.subplots(self.dim, self.dim)
            self.scatter_axes = np.reshape(self.scatter_axes, (self.dim, self.dim))
            for i in range(self.dim):
                for j in range(self.dim):
                    _j_margin = (self.upper_bound[j] - self.lower_bound[j]) * 0.1
                    _i_margin = (self.upper_bound[i] - self.lower_bound[i]) * 0.1
                    if j < i:
                        self.scatter_axes[i, j].set_xlim(
                            [self.lower_bound[j] - _j_margin, self.upper_bound[j] + _j_margin])
                        self.scatter_axes[i, j].set_ylim(
                            [self.lower_bound[i] - _i_margin, self.upper_bound[i] + _i_margin])
                        if i == self.dim - 1:
                            self.scatter_axes[i, j].set_xlabel(self.state_names[j])
                        if j == 0:
                            self.scatter_axes[i, j].set_ylabel(self.state_names[i])
                        self.scatter_axes[i, j].grid(True)
                        self.scatter_axes[i, j].set_aspect((self.upper_bound[j]
                                                            - self.lower_bound[j]
                                                            + 2 * _j_margin)
                                                           / (self.upper_bound[i]
                                                              - self.lower_bound[i]
                                                              + 2 * _i_margin))
                    elif i == j:
                        self.scatter_axes[i, j].set_xlim(
                            [self.lower_bound[j] - _j_margin, self.upper_bound[j] + _j_margin])
                        self.scatter_axes[i, j].set_ylim([-0.1, 1.1])
                        if i == self.dim - 1:
                            plt.xlabel(self.state_names[i])
                        self.scatter_axes[i, j].grid(True)
                        self.scatter_axes[i, j].set_aspect((self.upper_bound[j]
                                                            - self.lower_bound[j]
                                                            + 2 * _j_margin)
                                                           / (1.2))
                    elif j > i:
                        self.scatter_axes[i, j].remove()

        for i in range(self.dim):
            for j in range(self.dim):
                if j < i:
                    if hasattr(self, "_last_sample"):
                        self.scatter_axes[i, j].plot(self._last_sample[j], self._last_sample[i], ".", color="blue")
                    self.scatter_axes[i, j].plot(online_data[j], online_data[i], ".", color="orange")

                if j == i:
                    self.scatter_axes[i, j].plot(self.coverage_data[j], np.zeros_like(self.coverage_data[i]), ".",
                                                 color="blue")
                    self.scatter_axes[i, j].plot(self.suggested_sample[i], 0, ".", color="orange")


        self._last_sample = self.suggested_sample

        plt.pause(0.001)

    def plot_scatter(self,
                     scatter_kwargs={}):


        if not hasattr(self, "scatter_fig") or not self.render_online:
            self.scatter_fig, self.scatter_axes = plt.subplots(self.dim, self.dim)
            self.scatter_axes = np.reshape(self.scatter_axes, (self.dim, self.dim))

            for i in range(self.dim):
                for j in range(self.dim):
                    _j_margin = (self.upper_bound[j] - self.lower_bound[j]) * 0.1
                    _i_margin = (self.upper_bound[i] - self.lower_bound[i]) * 0.1

                    if j < i:
                        self.scatter_axes[i, j].set_xlim([self.lower_bound[j] - _j_margin,
                                                          self.upper_bound[j] + _j_margin])
                        self.scatter_axes[i, j].set_ylim([self.lower_bound[i] - _i_margin,
                                                          self.upper_bound[i] + _i_margin])
                        if i == self.dim - 1:
                            self.scatter_axes[i, j].set_xlabel(self.state_names[j])
                        if j == 0:
                            self.scatter_axes[i, j].set_ylabel(self.state_names[i])

                        self.scatter_axes[i, j].grid(True)
                        self.scatter_axes[i, j].set_aspect((self.upper_bound[j]
                                                            - self.lower_bound[j]
                                                            + 2 * _j_margin)
                                                           / (self.upper_bound[i]
                                                              - self.lower_bound[i]
                                                              + 2 * _i_margin))
                    elif i == j:
                        self.scatter_axes[i, j].set_xlim(
                            [self.lower_bound[j] - _j_margin, self.upper_bound[j] + _j_margin])
                        self.scatter_axes[i, j].set_ylim([-0.1, 1.1])
                        self.scatter_axes[i, j].grid(True)
                        self.scatter_axes[i, j].set_aspect((self.upper_bound[j]
                                                            - self.lower_bound[j]
                                                            + 2 * _j_margin)
                                                           / (1.2))

                        if i == self.dim - 1:
                            plt.xlabel(self.state_names[i])

                    elif j > i:
                        self.scatter_axes[i, j].remove()

        for i in range(self.dim):
            for j in range(self.dim):
                if j < i:
                    scatter_kwargs.setdefault("color", "blue")
                    scatter_kwargs.setdefault("marker", ".")
                    scatter_kwargs.setdefault("linestyle", "")
                    self.scatter_axes[i, j].plot(self.coverage_data[j], self.coverage_data[i], **scatter_kwargs)

                if j == i:
                    kde = kale.KDE(self.coverage_data[i], bandwidth=self.bandwidth)
                    points, density = kde.density(np.linspace(self.lower_bound[i], self.upper_bound[i], 500),
                                                  probability=True)
                    self.scatter_axes[i, j].plot(points, density)
                    self.scatter_axes[i, j].set_ylim([-0.1, np.max([1.1, np.max(density)])])
                    _j_margin = (self.upper_bound[j] - self.lower_bound[j]) * 0.1
                    self.scatter_axes[i, j].set_aspect((self.upper_bound[j]
                                                        - self.lower_bound[j]
                                                        + 2 * _j_margin)
                                                       / (np.max([1.1, np.max(density)]) - 0.1))
                    plt.grid(True)

        plt.show()

        return self.scatter_fig, self.scatter_axes
