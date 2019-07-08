using ProximalAlgorithms

const ForwardBackwardSolver = Union{
    ProximalAlgorithms.ForwardBackward,
    ProximalAlgorithms.ZeroFPR,
    ProximalAlgorithms.PANOC,
}

const default_solver = ProximalAlgorithms.PANOC
