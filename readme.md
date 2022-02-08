1. Lee to prepare ERSEM restarts from iMARNET (ORCA1 --> eORCA1) and then put WOA (N,P,DO) and GLODAP (DIC, TA)

2. LEE to prepare dust, CO2, N2O, Ndep from IMARNET and other sources

3. Giovanni to test  on full ARCHER2 and share with Lee (including physics restart and grid)

4. Giovanni to move from 3 tracers to full ERSEM

5. Yuri to ask James/Dale for ady

6. Yuri to share ERSEM restart with LEE and ARCHER2 link



Notes from meeting - 2021-11-25:
      Mission Atlantic - 3D simulation of global nemo-ersem
        NOC has delievered a basic 1 degree setup
        working on finalising that, moving to 1/.4 degree
        adding tides
        provenance
        they're doing development
        we should start running the physics set up so that we can see how it behaves globally
        Then see what changs we need to do, what kind of re-parameterise we need to do.
        
        Giovanni has put FABM into NEMO set up.
        Starting to set up NEMO-FABM.
        prepare NEMO-FabM-ERSEM for first simulation.
        today - plan and timeline for how to proceed.
        
        Giovanni - ARCHER run with NEMO-FABM. 1 month.
                passive tracers - a few months, but it works.
                soon upgrading to full ERSEM.
        
        progress on Hindcast:
                physics & forcing done.,
                rivers
                N deposition
                Fe deposition?
                ADY from iop: shading 
                sea-ice - not from the model? initial conditions
                atmposheric CO2/N2O from Mauna luoa, NOAA
                N20 not in Standard ERSEM.
                        could use atmospheric forcing N2O in water column       
                                feedback back to atmosphere?
                        air sea exchange of N2O.
        
        Hindcast:
                spin up?
                1 degree spin up with 
                ocean only 
                initial conditions: CMEMS + WOA/GLODAP (for carbon)
                try to start close to observations.
                
                starting year: 1980-2020 ish.
        
        Run for a year or something.
        
        ARCHER
        
        ARCHER2 has two different systems. 
                small cabinet system - slowand small.
                new bigger system very recently started.
                similar hardware but different software libraries.
                Full system is currently free.
        
        Following James Harles set up - which comes with scripts that automatically do lots of stuff.
        
        Giovanni - 
                has a set up ready to go - ish.
        
        Action & Jobs:
                Get set up on ARCHER2.
                Start moving ERSEM into the run.
                Start getting the forcing 


