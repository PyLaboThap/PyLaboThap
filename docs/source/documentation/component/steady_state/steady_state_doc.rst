Steady-State Models
===================

This section provides an overview of the steady-state models for components in the PyLaboThap library.

Steady-state models serve as the fundamental building blocks for simulating systems. They are used to represent the behavior of 
individual components under steady conditions and can be interconnected using connectors to model more complex systems. Components 
can be defined in two different ways:

Connector Approach
-------------------
In this approach, the inputs for components are set through connectors. This is particularly useful when multiple components are linked 
to form a larger system. Connectors facilitate the transfer of inputs and outputs between components, enabling the simulation of system-wide 
behavior.

.. image:: ../../../../figures/component/Component_connectors.png
   :alt: Components connectors description.
   :width: 100%


Input/Output Approach
---------------------
In this approach, inputs are set directly within each component. This is more suitable when using components individually, without linking 
them to other components.

.. image:: ../../../../figures/component/Component_in_out.png
   :alt: Components inlet/outlet description.
   :width: 100%

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   base_component/base_component_doc
   models/models_doc

