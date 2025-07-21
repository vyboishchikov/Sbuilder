<h1> Sbuilder &ndash; a Python code to successively attache benzene rings to a given input geometry (xyz coordinates) </h1>

<h2>Introduction and User's Guide</h2>

<p><i>Sbuilder</i> is a Python code that adds benzene rings to any available C&ndash;C bond of an aromatic system. 
The program looks for all available C&ndash;C bonds (effectively for H&ndash;C&ndash;C&ndash;H fragments), eliminates
the hydrogen atoms and fuses a new benzene ring to this C&ndash;C bond. The resulting structure is checked for 
uniqueness and, optionally, optimized using a force field.

<p> The <a href="https://github.com/vyboishchikov/Sbuilder/blob/main/SBuilder.py">Sbuilder Python code</a> can be downloaded.
To run the program from the command line, use the following syntax:</p>
<code>Sbuilder.py <i>xyz-file</i> -gen number_of_generations -out <i>Output.xyz</i> [-opt uff|mmff]</code>

<p><i><b>Warning:</b></i> The <i>xyz-file</i> should contain atomic symbols and Cartesian coordinates (in &#8491;). 
Do not include any header in the file.</p>

<p> The program needs the following Python modules: <code>Numpy</code>, <code>Random</code>, <code>NetworkX</code>.

<p>The code was used in the following publication:<br>
<p>I. Sarfraz, S. F. Vyboishchikov, M. Sol&agrave;, A. Artigas, <i>submitted</i>, <b>2025</b>.

<p>For questions related to the <i>Sbuilder</i> program, please contact
<a href="mailto:vyboishchikov@googlemail.com">Sergei Vyboishchikov</a>.</p>
</body>
</html>
