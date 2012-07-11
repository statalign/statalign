package com.ppfold.main;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class AlignmentReader {

	public static Alignment readAlignment(String fullFileName) throws Exception{

		try {
			FileReader input = new FileReader(fullFileName);
			BufferedReader bufRead = new BufferedReader(input);
			String line;
			List<String> lines = new ArrayList<String>(); 

			line = bufRead.readLine();
			if(line!=null){
				lines.add(line.toString());
			}

			// Read through file one line at time. 
			while (line != null){
				line = bufRead.readLine();
				if(line!=null){
					lines.add(line.toString());
					//System.out.println(line);
				}
			}
			bufRead.close();

			Alignment align = parse(lines);

			return align;
		}
		catch (IOException e){
			throw new Exception("Error reading alignment file! Check that the name is OK.");
			
		}
		catch(Exception e){
			throw new Exception("Error parsing the alignment file!");
		}
	}
	
	public static Alignment readAlignmentFromStringList(List<String> lines) throws Exception{
		Alignment align = parse(lines);
		return align;
	}
	
	
	public static Alignment parse(List<String> lines) throws Exception{
		
		try{
			List<String> sequences = new ArrayList<String>();
			List<String> names = new ArrayList<String>();

			int i = 0; //line counter

			String line = "";
			String name = "";
			String sequence = ""; 
			while(i<lines.size()){
				line = lines.get(i);
				if(line.startsWith(">")){
					//if the line starts with > then we have found a sequence 
					//read name
					name = line.substring(1);
					i++;
					//fill the sequence corresponding to this name
					while(i<lines.size()){
						line = lines.get(i);				
						if(!line.startsWith(">")){
							line = line.replace(" ", "");
							sequence = sequence.concat(line.toLowerCase().trim());
							i++;
						}
						else{
							i--;
							break;
						}
					}
					names.add(name.trim());
					//	System.out.println("Added sequence name: " + name.trim());
					sequences.add(sequence);
					//	System.out.println("Added sequence: " + sequence);
					name="";
					sequence="";

				}
				i++;
			}

			return new Alignment(sequences, names);
		}
		catch(Exception e){
			throw new Exception("Error: The alignment could not be parsed. Check input.");
		}
	}

	
}

