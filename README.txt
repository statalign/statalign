VARNA is a tool for the automated drawing, visualization and annotation of the 
secondary structure of RNA, designed as a companion software for web servers 
and databases.
Copyright (C) 2008  Kevin Darty, Alain Denise and Yann Ponty.
electronic mail : Yann.Ponty@lri.fr
paper mail : LRI, bat 490 Université Paris-Sud 91405 Orsay Cedex France

The latest version of this software can be found at:
  http://varna.lri.fr

%%%%%%%%%%%%%%%%%%%%%%%%% INVOKING THE GUI VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A basic editor demonstrating the features of VARNA can be spawned through the 
following command:

java -jar VARNAv3-1.jar applis.VARNAGUI 

%%%%%%%%%%%%%%%%%%%%%% INVOKING COMMAND LINE VERSION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

VARNA can also be used as a standalone application for producing pictures. For 
such an application, the VARNAcmd class was coded, and should be invoked in 
the following way:

   java -jar VARNAv3-1.jar applis.VARNAcmd 
     [-i inputFile|-sequenceDBN XXX -structureDBN YYY] -o outFile [opts]

Where:
    * inFile: An input file using one of the supported formats (Vienna, CT, BPSeq or RNAML).
    * XXX: An RNA sequence.
    * YYY: A dot-bracket representation of this RNA.
    * outFile: An output file whose format is guessed from the extension.

Many additional options are supported, matching the syntax of the Applet's parameters. 
Namely, if param is a valid parameter for the VARNA Applet, then -param is an option
that has the exact same arguments and effects in the standalone, command-line version.

%%%%%%%%%%%%%%%%%%%%%% LEGAL STUFF AND USUAL DISCLAIMER %%%%%%%%%%%%%%%%%%%%%%%%%%

 This file is part of VARNA version 3.1.
 VARNA version 3.1 is free software: you can redistribute it and/or modify it 
 under the terms of the GNU General Public License as published by the Free 
 Software Foundation, either version 3 of the License, or (at your option) any 
 later version.

 VARNA version 3.1 is distributed in the hope that it will be useful, but WITHOUT 
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
 A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with 
 VARNA version 3.1 (See LICENSE.txt file). If not, see http://www.gnu.org/licenses.
