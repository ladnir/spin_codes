"""Extra proposal search for pairwise routing; outward verification unchanged."""
from fractions import Fraction as Q
import scalar_cover
import regional_count


def retune(model,cell,base,factors=(Q(3,4),Q(1,2),Q(5,4),Q(3,2),Q(7,8),Q(9,8))):
    score,witness=base
    stop=scalar_cover.proposal_cutoff(model)
    if score<stop or 'regional_count_parts' not in witness:return base
    lam=Q(witness['parameters'][0]);best=base
    for factor in factors:
        if Q(factor)<=0:raise ValueError('positive output-tilt factor required')
        trial=dict(witness,parameters=[str(lam*Q(factor)),*witness['parameters'][1:]])
        candidate=regional_count.propose(model,cell,trial)
        print('PAIRWISE REGIONAL RETUNE factor',str(factor),'proposal',candidate[0],flush=True)
        if candidate[0]<best[0]:best=candidate
        if best[0]<stop:break
    return best


class Model(scalar_cover.Model):
    def proposal(self,cell):
        return retune(self,cell,super().proposal(cell))
