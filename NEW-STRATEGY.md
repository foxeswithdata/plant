# Documenting 'plant' strategy extension

This document is intended for people wanting to extend the 'plant' package with a new strategy.
It is meant for people familiar with coding and reading code but perhaps new to C++ or needing a refresher.
I will cover some basic C++ content with the intention of naming the C++ idiosynchracies so that the reader can further research those.

The R-C++ interface is not covered in this document as that interface should be generated automatically when the package is compiled.
The R usage of the package is not covered either and the user should refer to the documentation provided by Falster et al 2016.

In order to compile the package use the Makefile that comes with the original package.
There should be no need to edit the contents of the Makefile.

New strategy will be referred to throughout the code and this document as 'NEW_strategy'.
When working with the code change all those references to your preferred strategy name, including the file names. e.g. NEW_strategy -> ES20_strategy.

Below is the list of files that might need to be changed:

* `inst/include/plant.h`
* `inst/include/plant/model/NEW_strategy.h`
* `inst/include/plant/model/NEW_environment.h` (if changing environment)
* `inst/include/plant/model/Assimilator.h` (if changing photosynthesis. not covered here. yet?)
* `src/NEW_strategy.cpp`


The h files are header files and will include the definitions of classes, methods and parameters.
The cpp files contain the implementation of those methods.
In some shorter and simple cases the header file might have some simple or important methods but mostly it's important to implement them in cpp files.

The list is not exclusive and depending on user needs there might be requirement to further change things.
A good way to do that might be to extend a class with your own and use that one in the model implementation.

The way I organise this document is through individual files.
It follows the logic of setting up the changes to be part of the project, followed by implementing the main elements of the new strategy which will be contained in the NEW_strategy files and finally I make a few notes about the environment.
I cover these file by file in terms of the implementation.

## inst/include/plant.h

This is the header class for the model definition. Contains links to all the header files in the model.

Add `#include <filepath/filename.h>` for any new C++ class or file.

There are commented out links for `NEW_strategy.h` and `NEW_environment.h`.
These can be uncommented and file names can be changed and/or other files can be added as becomes needed.

## 'inst/include/plant/model/NEW_strategy.h'

The NEW_strategy.h is a template header file that can be filled out with model and strategy details.
The template code provides clues and strategies and some explanations about what individual methods do in the context of the package.
It's also useful to look at the existing FF16_Strategy file to get an idea of some details.

This class inherits from the Strategy.h class so it's important to implement all methods from that class.
I don't cover all of those methods in this document so please refer to Strategy.h and FF16_Strategy.h when editing.

Here are some steps that can be followed in the initial instance:

1. Replace definition at top with strategy name (lines 2 and 3) (this has to be done for all the new header files)
2. Make sure to include all necessary dependencies (including the environment used; replace NEW_environment.h (line 11) with new or existing strategy).
Might need to include other dependencies as need comes up.
3. Change class definition
  `class NEW_Strategy: public Strategy<NEW_Environment>`
  Should include template of `Strategy<Environment_used>`.
  Potentially could inherit existing strategy eg.
  `class NEW_Strategy: public FF16_Strategy, Strategy<NEW_Environment>`
  but I still need to figure out how that would work further.
  If you're unsure of this set-up look up templates and inheritance in C++ to further understand these details.
4. *States* are represented as a vector rather than as individual variables. Define names in `state_names()` method.
State values are themselves kept in `Internals` object `vars` which is passed to the relevant methods (not stored in the strategy file).
In the strategy file need to store both the names and indices of strategies in a map file (`state_index`).
Need to also make sure that `state_size()` returns correct number of states.
5. *Auxiilary states* work in the same way as described above for states.
(TODO: what is the purpose of auxilary states?)
6. The method `compute_rates()` contains all the necessary calls of behaviour of the strategy also all incoming information (eg. in the environment object).
The method is inherited from the Strategy class and so the definition shouldn't be changed for the model to work.
7. Need to add any helper methods. Preferably one method for each equation.

Note on assimilation. This is implemented in the Assimilation.h class.
This modularises the photosynthesis behaviour and any photosynthesis changes should be implemented there.

## 'src/NEW_strategy.cpp'

In this file the implementation of the

This is the implementation of the strategy file as described above.
Again, most notes are included in the code.

Here I just point out some important bits.

1. Most important here is the `compute_rates()` method.
Recover all the necessary states through the Internals vars: `vars.state(STATE_1_INDEX)` or
`vars.state(state_index.at("state_name"))`.
Recover auxilary states in a similar fashion using `vars.aux()`.
To set the rates of change of states use `vars.set_rate(index_num, state_dt);`.
Set auxilary states through `vars.set_aux(aux_index.at("aux_name"), aux_value)`.
2. Implement any of the declared model helper methods.
Also implement all the methods from the Strategy class.
3. In `prepare_strategy()` do any initial extra setups for the strategy such as initialising the control and the assimilator.
Control can be disregarded for now as it contains some non-biological parametrisation (more info on that in the Control.h class).

## 'inst/include/plant/model/NEW_environment.h'

TODO: compile the `NEW_environment.h` template class.
Any new environment changes will be implemented here.
There's a chance that there will be no need to edit the environment but if so then here is some guidelines.

The new environment class will inherit from `Environment` class.
It will hold information on the disturbance, time, light environment etc.
There's potential to also extend from another environment class like `FF16_Environment`.

1. The environment holds a height-distributed variation of the light environment.
This is held in the adaptive_interpolator object called `environment_interpolator`.
A function is passed when computed (in the c`ompute_environment()` method) which which calculates the intermediate steps between bounds.
This behaviour is held in the `AdaptiveInterpolator.h` class.
The actual function that computes this is actually defined elsewhere and passed to the environment so the environment is really only a place-holder for this (except for the `K_I` light decay parameter).
NOTE: this is mostly based on the previous iteration and will need to be probably
updated as the code has changed.
2. There are multiple `template <typename Function>` methods in the original and new Environment class.
These allow for a function to be passed and used in the method.
If you want to find out what is normally passed through to here you will need to find what function is passed into an environment object.
3. TODO: I ignore the disturbance regime for now as I don't use it. I need to figure out how to switch it off. TODO
4. TODO: figure out and discribe the seed regime.
5. TODO: Important question. How is the time actually changed?

## Further Work

This is a living document that I will edit as I proceed with answering my own questions and adding solutions.
This document is an introduction into the process and will be updated with individual problems and how they can be changed, thereby changing the format of this document.
