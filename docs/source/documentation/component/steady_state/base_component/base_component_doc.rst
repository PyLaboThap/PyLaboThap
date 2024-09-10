BaseComponent Class
-------------------
The `BaseComponent` class is the parent class for all components in the PyLaboThap library. All components inherit the properties and methods of the `BaseComponent` class. Key methods provided by this class include:

.. autoclass:: component.base_component.BaseComponent

Example of a Steady-State Component Model
------------------------------------------

The following example demonstrates how to create a compressor model using the PyLaboThap library. This example is based on the model (XXX). For more details on this model, refer to the documentation [here](lien).

Inputs for the model can be defined using either the connector approach or the input/output approach. The figures below illustrate these two approaches for the compressor model:

.. raw:: html

   <div class="side-by-side">

.. image:: ../../../../../figures/component/compressor_connectors.png
   :alt: Connectors approach for a compressor model.
   :width: 50%

.. image:: ../../../../../figures/component/compressor_in_out.png
   :alt: Input/Output approach for a compressor model.
   :width: 50%

.. raw:: html

   </div>

For each model, the following methods need to be implemented:

- **Connectors**: Specify connectors in the `__init__` method of the class.
  
  .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
     :language: python
     :lines: 22-27

- **Required Inputs**: Define the required inputs for the model to work in the `get_required_inputs` method.
  
  .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
     :language: python
     :lines: 29-32

- **Sync Inputs**: In the `sync_inputs` method, the inputs dictionary (containing all of the inputs) is synchronized with the connectors' states. This method searches for inputs inside the connectors if inputs are set through them.
  
  .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
     :language: python
     :lines: 34-47

- **Set Inputs**: In the `set_inputs` method, the inputs are set directly by the user. The connectors will thus be filled by the inputs.
  
  .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
     :language: python
     :lines: 49-65

- **Get Required Parameters**: List all the parameters required to run the model in the `get_required_parameters` method.
  
  .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
     :language: python
     :lines: 67-70

- **Print Setup**: Print the different connectors and inputs needed to run the simulation in the `print_setup` method.
  
  .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
     :language: python
     :lines: 67-70

- **Solve**: Implement the model-solving logic in the `solve` method.

- **Update Connectors**: In the `update_connectors` method, update the connectors after calculation.
  
  .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
     :language: python
     :lines: 339-348

- **Print Results**: Use the `print_results` method to print the results calculated by the model.
  
  .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py
     :language: python
     :lines: 350-357

- **Print States of Connectors**: Print the connectors' state in the `print_states_connectors` method.
  
  .. literalinclude:: ../../../../../../library/component/steady_state/volumetric_machine/compressor/semi_empirical/simulation_model.py


More information on this model can be found [here](lien).


