import numpy as np
from scipy import stats
import ddm

class OverlayNonDecisionGaussian(ddm.Overlay):
    """ Courtesy of the pyddm cookbook """
    name = "Add a Gaussian-distributed non-decision time"
    required_parameters = ["ndt_location", "ndt_scale"]

    def apply(self, solution):
        # Extract components of the solution object for convenience
        corr = solution.corr
        err = solution.err
        dt = solution.model.dt
        # Create the weights for different timepoints
        times = np.asarray(list(range(-len(corr), len(corr)))) * dt
        weights = stats.norm(scale=self.ndt_scale, loc=self.ndt_location).pdf(times)
        if np.sum(weights) > 0:
            weights /= np.sum(weights)  # Ensure it integrates to 1
        newcorr = np.convolve(weights, corr, mode="full")[len(corr):(2 * len(corr))]
        newerr = np.convolve(weights, err, mode="full")[len(corr):(2 * len(corr))]
        return ddm.Solution(newcorr, newerr, solution.model,
                            solution.conditions, solution.undec)


class DriftTtaDistanceStatic(ddm.models.Drift):
    name = "Drift depends on initial TTA and distance but doesn't change over time"
    required_parameters = ["alpha", "beta", "theta"]
    required_conditions = ["tta_condition", "d_condition"]

    def get_drift(self, t, conditions, **kwargs):
        return self.alpha * (conditions["tta_condition"] + self.beta * conditions["d_condition"] - self.theta)


class DriftTtaDistanceDynamic(DriftTtaDistanceStatic):
    name = "Drift dynamically depends on the real-time values of TTA and distance"

    def get_drift(self, t, conditions, **kwargs):
        return self.alpha * (conditions["tta_condition"] - t + self.beta
                             * (conditions["d_condition"] - t * conditions["d_condition"] / conditions["tta_condition"])
                             - self.theta)


class BoundCollapsingTta(ddm.models.Bound):
    name = "Bounds dynamically collapsing with TTA"
    required_parameters = ["b_0", "k", "tta_crit"]
    required_conditions = ["tta_condition", "d_condition"]

    def get_bound(self, t, conditions, **kwargs):
        tau = conditions["tta_condition"] - t
        return self.b_0 / (1 + np.exp(-self.k * (tau - self.tta_crit)))

class BoundCollapsingTtaDistance(ddm.models.Bound):
    name = "Bounds dynamically collapsing with TTA and distance"
    required_parameters = ["b_0", "k", "beta_boundary", "theta_boundary"]
    required_conditions = ["tta_condition", "d_condition"]

    def get_bound(self, t, conditions, **kwargs):
        tta = conditions["tta_condition"] - t
        d = (conditions["d_condition"] - t * conditions["d_condition"] / conditions["tta_condition"])
        return self.b_0 / (1 + np.exp(-self.k * (tta + self.beta_boundary*d - self.theta_boundary)))

class BoundCollapsingGeneralizedGap(ddm.models.Bound):
    name = "Bounds dynamically collapsing with generalized gap"
    required_parameters = ["b_0", "k", "beta", "theta"]
    required_conditions = ["tta_condition", "d_condition"]

    def get_bound(self, t, conditions, **kwargs):
        tta = conditions["tta_condition"] - t
        d = (conditions["d_condition"] - t * conditions["d_condition"] / conditions["tta_condition"])
        generalized_gap = (tta + self.beta * d)

        return self.b_0 / (1 + np.exp(-self.k * (generalized_gap - self.theta)))

class ModelStaticDriftFixedBounds():
    # simplest model, vanilla DDM with drift dependent on TTA and distance
    T_dur = 7.0
    # T_dur = 4.0 USED for training output_all_RT
    # T_dur = 3.0
    param_names = ["alpha", "beta", "theta", "b", "ndt_location", "ndt_scale"]

    def __init__(self):
        self.overlay = OverlayNonDecisionGaussian(ndt_location=ddm.Fittable(minval=0, maxval=2.5),
                                                  ndt_scale=ddm.Fittable(minval=0.001, maxval=1.0))
        self.drift = DriftTtaDistanceStatic(alpha=ddm.Fittable(minval=0.1, maxval=5),
                                            beta=ddm.Fittable(minval=0, maxval=1.0),
                                            theta=ddm.Fittable(minval=4, maxval=60))
        self.bound = ddm.BoundConstant(B=ddm.Fittable(minval=0.1, maxval=5))

        self.model = ddm.Model(name="Static drift defined by initial TTA and d, constant bounds",
                               drift=self.drift, noise=ddm.NoiseConstant(noise=1), bound=self.bound,
                               overlay=self.overlay, T_dur=self.T_dur)


class ModelDynamicDriftFixedBounds(ModelStaticDriftFixedBounds):
    def __init__(self):
        super().__init__()

        self.drift = DriftTtaDistanceDynamic(alpha=ddm.Fittable(minval=0.1, maxval=5.0),
                                             beta=ddm.Fittable(minval=0, maxval=1.0),
                                             theta=ddm.Fittable(minval=4, maxval=60))

        self.model = ddm.Model(name="Dynamic drift defined by real-time TTA and d, constant bounds",
                               drift=self.drift, noise=ddm.NoiseConstant(noise=1), bound=self.bound,
                               overlay=self.overlay, T_dur=self.T_dur)


class ModelDynamicDriftCollapsingBoundsTta(ModelDynamicDriftFixedBounds):
    param_names = ["alpha", "beta", "theta", "b_0", "k", "tta_crit", "ndt_location", "ndt_scale"]

    def __init__(self):
        super().__init__()

        self.bound = BoundCollapsingTta(b_0=ddm.Fittable(minval=0.5, maxval=5),
                                        k=ddm.Fittable(minval=0.1, maxval=2),
                                        tta_crit=ddm.Fittable(minval=3, maxval=6))

        self.model = ddm.Model(name="Dynamic drift defined by real-time TTA and d, bounds collapsing with TTA",
                               drift=self.drift, noise=ddm.NoiseConstant(noise=1), bound=self.bound,
                               overlay=self.overlay, T_dur=self.T_dur)

class ModelDynamicDriftCollapsingBoundsTtaDistance(ModelDynamicDriftCollapsingBoundsTta):
    param_names = ["alpha", "beta", "theta", "b_0", "k", "beta_boundary", "theta_boundary", "ndt_location", "ndt_scale"]

    def __init__(self):
        super().__init__()

        self.drift = DriftTtaDistanceDynamic(alpha=ddm.Fittable(minval=0.1, maxval=5.0),
                                             beta=ddm.Fittable(minval=0, maxval=1.0),
                                             theta=ddm.Fittable(minval=4, maxval=60))

        self.bound = BoundCollapsingTtaDistance(b_0=ddm.Fittable(minval=0.5, maxval=5),
                                                k=ddm.Fittable(minval=0.1, maxval=2),
                                                beta_boundary=ddm.Fittable(minval=0, maxval=1.0),
                                                theta_boundary=ddm.Fittable(minval=4, maxval=60))

        self.model = ddm.Model(name="Dynamic drift defined by real-time TTA and d, bounds collapsing with TTA and distance",
                               drift=self.drift, noise=ddm.NoiseConstant(noise=1), bound=self.bound,
                               overlay=self.overlay, T_dur=self.T_dur)

class ModelDynamicDriftCollapsingBoundsGeneralizedGap(ModelDynamicDriftCollapsingBoundsTta):
    param_names = ["alpha", "beta", "theta", "b_0", "k", "ndt_location", "ndt_scale"]

    def __init__(self):
        super().__init__()

        beta = ddm.Fittable(minval=0, maxval=1.0)
        theta = ddm.Fittable(minval=4, maxval=60)

        self.drift = DriftTtaDistanceDynamic(alpha=ddm.Fittable(minval=0.1, maxval=5.0),
                                             beta=beta,
                                             theta=theta)

        self.bound = BoundCollapsingGeneralizedGap(b_0=ddm.Fittable(minval=0.5, maxval=5),
                                                k=ddm.Fittable(minval=-0.5, maxval=2),
                                                beta=beta,
                                                theta=theta)

        self.model = ddm.Model(name="Dynamic drift defined by real-time TTA and d, bounds collapsing with generalized gap",
                               drift=self.drift, noise=ddm.NoiseConstant(noise=1), bound=self.bound,
                               overlay=self.overlay, T_dur=self.T_dur)

