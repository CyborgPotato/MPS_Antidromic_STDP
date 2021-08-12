

# Implementation of Antidromic Motorunit Model

This is a basic implementation of a compartmental antidromic motor
unit model based off the work of
<sup id="2f569add2af6f09765def7d7d132f335"><a href="#cisi08_simul_system_spinal_cord_motor" title="Rogerio Cisi \&amp; Andr\'e Kohn, Simulation System of Spinal Cord Motor Nuclei and  Associated Nerves and Muscles, in a Web-Based  Architecture, {Journal of Computational Neuroscience}, v(3), 520-542 (2008).">cisi08_simul_system_spinal_cord_motor</a></sup> modified to include STDP
for long term plasticity based off the work of
<sup id="e59e69cd72dd9b77907b4df914aa773e"><a href="#pedrosa17_role_neurom_cortic_plast" title="Victor Pedrosa \&amp; Claudia Clopath, The Role of Neuromodulators in Cortical  Plasticity. a Computational Perspective, {Frontiers in Synaptic Neuroscience}, v(nil), nil (2017).">pedrosa17_role_neurom_cortic_plast</a></sup>.


## Running:

To use the model refer to `main.py`.

**Important:** Neurons are connected in reverse
i.e. `neuron1.connect(neuron2)` connects `neuron1` and `neuron2` such
that `neuron2` synapses onto `neuron1`'s dendrite.

