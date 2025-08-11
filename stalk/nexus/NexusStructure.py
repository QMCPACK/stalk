#!/usr/bin/env python3

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from structure import Structure
from simulation import Simulation

from stalk.params.ParameterStructure import ParameterStructure


class NexusStructure(ParameterStructure):
    _jobs: list[Simulation] = None
    _sigma = None

    @property
    def jobs(self):
        return self._jobs
    # end def

    @jobs.setter
    def jobs(self, jobs):
        if jobs is None or len(jobs) == 0:
            self._jobs = None
        else:
            for job in jobs:
                if not isinstance(job, Simulation):
                    raise TypeError("Nexus job must be inherited from Simulation class!")
                # end if
            # end for
            self._jobs = jobs
        # end if
    # end def

    @property
    def generated(self):
        return self.jobs is not None
    # end def

    @property
    def finished(self):
        return self.generated and all(job.finished for job in self._jobs)
    # end def

    @property
    def analyzed(self):
        return self.finished and self.value is not None
    # end def

    @property
    def sigma(self):
        return self._sigma
    # end def

    @sigma.setter
    def sigma(self, sigma):
        if sigma is None:
            self._sigma = 1e100
        elif sigma >= 0:
            self._sigma = sigma
        else:
            raise ValueError(f'Sigma must be >= 0, provided: {sigma}')
        # end if
    # end def

    def get_nexus_structure(
        self,
        kshift=(0, 0, 0),
        **kwargs
    ):
        kwargs.update({
            'elem': self.elem,
            'pos': self.pos,
            'units': self.units,
        })
        if self.axes is not None:
            kwargs.update({
                'axes': self.axes,
                'kshift': kshift,
            })
        # end if
        return Structure(**kwargs)
    # end def

    def reset_value(self):
        super().reset_value()
        # Reset jobs upon value change
        self._jobs = None
    # end def

    def copy(
        self,
        **kwargs
        # params=None, params_err=None, label=None, pos=None, axes=None, offset=None
    ):
        tmp_jobs = self._jobs
        # Put jobs lists aside during copy
        self._jobs = None
        result = super().copy(**kwargs)
        # Recover jobs lists
        self._jobs = tmp_jobs
        return result
    # end def

# end class
