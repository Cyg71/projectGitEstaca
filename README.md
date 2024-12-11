# projectGitEstaca
Repository of the ESTACA project about GitHub. Developped by Valentine LEGAT and Cyriaque GUILLOT, ESTACA students.

## Introduction
This project is a J2 Secular Orbital propagator. It will calculate the successive position and velocity vectors in ECI, taking into account the gravity of the Earth (2-body perturbation) and the influence of the 1st of the oblateness of the Earth (J2 parameter). It computes the evolution of the orbital elements, especially the precession of the longitude of the Right Ascension of the Ascending Node (RAAN or Ω) and the precession of the Argument of perigee (ω).

## Background
The Earth is not a perfect sphere. Due to the centrifugal forces, the equatorial radius of the Earth is slightly bigger than its polar radius. This ununiformity causes a perturbation in the gravitational field that the 2-body solicitation does not take into account. To mitigate this deviation, the J2 parameter can be implemented in the propagation code. It is one in many coefficients that account for the oblateness of the Earth, and the most significant.

## Parameters
- TLE: In the TLE tuple, copy and paste the TLE of the orbital element of your choosing. The Two-Line Element sets are available on various website such as https://www.n2yo.com/
- nb_days: The number of days along which you want your simulation to propagate your orbit.
- dt: The time step of your simulation

## Features
- Propagates the satellite around the Earth, accounting for 2-body and J2 perturbations
- Plots the propagated trajectory in 3D around the Earth
- Plots the ECI coordinates X, Y and Z versus time

## Installation for reuse
1. To setup the environnement, make sure to install miniconda https://docs.anaconda.com/miniconda/miniconda-install/
2. Create your python environnement
    ```bash
    conda create -y -n my_env python=3.11
    conda activate my_env
    pip install pdm ruff==0.3.3 pre-commit
3. Go to your folder and install the packages using pdm. Initialize the project if needed
    ```bash
    cd ~/your_work_folder
    pdm init
    pdm install
    pip install -r requirements.txt
    cd src\projectgitestaca
    python orbitalPropagation.py

## Contact
Valentine LEGAT - Valentine.LEGAT@estaca.eu - ESTACA Paris-Saclay <br />
Cyriaque GUILLOT - Cyriaque.GUILLOT@estaca.eu - ESTACA Paris-Saclay

## Acknowledgements
This project was made possible thanks to the introduction and help of Nicolas Dagoneau, Researcher in Astrophysics and software engineer at CEA. <br />
We would also like to thank Arnaud Somville whose template was used to create and document our project. His help has been paramount in the development of our code. <br />
Template available here: https://github.com/arnaudsomville/python_project_template/tree/main