+++
date = "2021-06-17T11:00:13+01:00"
title = "OGS Conventions"
author = "Joerg Buchwald"
weight = 1

aliases = [ "/docs/userguide/",
            "/docs/quickstart/",
            "/docs/quickstart/basics/quickstart" ] # Redirect for Hydrology III tutorial

[menu.docs]
name = "User Guide"
identifier = "userguide"
weight = 1
post = "Download, install and run an OGS benchmark in 5 minutes! No development setup required."

[menu]
  [menu.userguide]
    parent = "basics"
+++

## Units

In general OpenGeoSys is not using any intrinsical units, i.e. all equations are [non-dimensionalized](https://en.wikipedia.org/wiki/Nondimensionalization).
However, there are some exemptions to the rule as for instance some empirical laws
like some expressions for the vapor diffusion that assume the temperature to be given in Kelvin.

List of exceptions:
 - ...

Dimensional analysis of source and boundary terms...

## Stress sign

OpenGeoSys uses the positive sign convention for tensile stresses, whereas compressional stresses carry a negative sign.

## Thermal expansivity

In OpenGeoSys the term themal expansivity corresponds to the volumentric property of a fluid and to the
linear thermal expansion coefficient if it is used in the context of a solid where anisotrpic expansion can also occur.

## Order of process variables/global components

OpenGeoSys uses internally a vector of primary vector components ("global component vector").
Although, which set of components is used, varies from process to process, for monolithic process schemes,
the order of primary variable components is maintained: $T$,$p$,$u_{x}$,$u_{y}$,$u_{z}$.
This order is used e.g., in the definition and the output of the per component convergence of the non-linear solver.

## Kelvin mapping

To map the elastic tensor OpenGeoSys internally uses a Kelvin mapping with and an adapted component ordering for computational reasons.
For 2D the Kelvin-Vector of the stress tensor looks like $\sigma=(\sigma_{xx},\sigma_{yy},\sigma_{zz},\sqrt(2)\sigma_{xy})$ whereas the 3D version reads as $\sigma=(\sigma_{xx},\sigma_{yy},\sigma_{zz},\sqrt(2)\sigma_{xy}, \sqrt(2)\sigma_{yz},\sqrt(2)\sigma_{xz})$.
