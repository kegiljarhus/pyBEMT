
class Fluid:
    def __init__(self, cfg):
        self.rho = cfg.getfloat('fluid','rho')
        self.mu = cfg.getfloat('fluid','mu')


