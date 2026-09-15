#!/usr/bin/env python3
'''TransitionStateSearch class for finding transition pathways.'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from pathlib import Path

from numpy import sort, dot, array

from stalk.io.stalk_path import StalkPath
from stalk.lsi.pathway_image import PathwayImage
from stalk.params.parameter_set import ParameterSet


class TransitionPathway(StalkPath):
    _images: list[PathwayImage] = []  # list of LineSearchIteration objects

    def __init__(
        self,
        images: list[ParameterSet] = None,
        path: str | Path | None = None,
    ):
        StalkPath.__init__(self, path=path)
        self._images = []
        if images is not None:
            # add image A
            self.add_image(images[0])
            # add image B
            self.add_image(images[-1])
            for image in images[1:-1]:
                self.add_image(image)
            # end for
        # end def
    # end def

    # Return a list of all pathway images
    @property
    def images(self):
        return self._images
    # end def

    # Return a list of intermediate pathway images
    @property
    def intermediate_images(self):
        return self._images[1:-1]
    # end def

    @property
    def pointA(self):
        if len(self) >= 1:
            return self.images[0]
        # end if
    # end def

    @property
    def pointB(self):
        if len(self) >= 2:
            return self.images[-1]
        # end if
    # end def

    @property
    def difference(self):
        if self.pointB is not None:
            return self.pointB.structure.params - self.pointA.structure.params
        # end if
    # end def

    @property
    def pathway_init(self):
        params = []
        params_err = []
        for image in self.images:
            params.append(image.structure_init.params)
            params_err.append(image.structure_init.params_err)
        # end for
        return array(params), array(params_err)
    # end def

    @property
    def pathway_final(self):
        params = []
        params_err = []
        for image in self.images:
            params.append(image.structure_final.params)
            params_err.append(image.structure_final.params_err)
        # end for
        return array(params), array(params_err)
    # end def

    def add_image(self, image: ParameterSet, rc=None):
        if self.pointA is None:
            # add point A
            pw_image = PathwayImage(
                image,
                reaction_coordinate=0.0,
                path=self.path / 'image_A'
            )
            self.images.append(pw_image)
        elif self.pointB is None:
            # add point B
            pw_image = PathwayImage(
                image,
                reaction_coordinate=1.0,
                path=self.path / 'image_B'
            )
            self.images.append(pw_image)
        else:
            if rc is None:
                rc = self._calculate_rc(image)
            # end if
            if rc <= 0.0:
                raise ValueError("Cannot add intermediate image with reaction coordinate <= 0")
            elif rc >= 1.0:
                raise ValueError("Cannot add intermediate image with reaction coordinate >= 1")
            else:
                # Insert next to last, presuming ordering by reaction coordinate
                pw_image = PathwayImage(
                    image,
                    reaction_coordinate=1.0,
                    path=self.path / f'image_{rc:+5.4f}/'
                )
                self.images.insert(-1, pw_image)
                sort(self.intermediate_images)
            # end if
        # end if
    # end def

    def calculate_hessians(
        self,
        **hessian_args,
    ):
        self.pointA.calculate_hessian(
            tangent=None,
            **hessian_args
        )
        self.pointB.calculate_hessian(
            tangent=None,
            **hessian_args
        )
        for i, image in enumerate(self.intermediate_images):
            im_prev = self.images[i]
            im_next = self.images[i + 2]
            tangent = im_next.structure.params - im_prev.structure.params
            image.calculate_hessian(
                tangent=tangent,
                **hessian_args
            )
        # end for
    # end def

    def generate_surrogates(self, **surrogate_args):
        for _, image in enumerate(self.images):
            image.generate_surrogate(**surrogate_args)
        # end for
    # end def

    def optimize_surrogates(self, **optimize_args):
        for _, image in enumerate(self.images):
            image.optimize_surrogate(**optimize_args)
        # end for
    # end def

    def run_linesearches(self, **lsi_args):
        for _, image in enumerate(self.images):
            image.run_linesearch(**lsi_args)
        # end for
    # end def

    def _calculate_rc(self, image: ParameterSet):
        rc = (dot(self.difference, image.params - self.pointA.structure.params) /
              dot(self.difference, self.difference))
        return rc
    # end def

    def __len__(self):
        return len(self.images)
    # end def

# end class
