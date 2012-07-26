package statalign.distance;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public 

class AlignmentSample {

	final static String pathToSamples = "/home/ingolfur/Dropbox/RNA and StatAlign/TestRNAData/mpd";

	public static void main(String[] args){
		File reffolder = new File(pathToSamples);
		File[] listOfRefFiles = reffolder.listFiles();
		ArrayList<AlignmentSample> yeah = new ArrayList<AlignmentSample>();
		for(File file  : listOfRefFiles){
			AlignmentSample temp = AlignmentSample.loadAlignments(file);
			 yeah.add(temp);
			 for(int i = 0; i<temp.posteriors.size(); i++){
				 System.out.println(temp.posteriors.get(i)); 
			 }
			 
			 System.out.println(file.getName());
		}
		
	}

	Alignment reference = new Alignment();
	Alignment mpd = new Alignment();
	ArrayList<Double> posteriors = new ArrayList<Double>();
	ArrayList<Alignment> samples = new ArrayList<Alignment>();

	public String toString() {
		String ret = "";
		ret += "%reference\n" + reference.toString();
		for (int i = 0; i < samples.size(); i++) {
			ret += "%" + i + "\n";
			ret += samples.get(i);
		}
		ret += "%mpd\n" + mpd.toString();
		ret += "%posteriors\n";
		for (int i = 0; i < posteriors.size(); i++) {
			ret += posteriors.get(i) + "\n";
		}
		return ret;
	}


	public static AlignmentSample loadAlignments(File file) {

		AlignmentSample sample = new AlignmentSample();
		try {
			BufferedReader buffer = new BufferedReader(new FileReader(file));
			String textline = null;

			String alignmentString = "";
			String name = null;

			ArrayList<String> sequences = new ArrayList<String>();
			ArrayList<String> sequenceNames = new ArrayList<String>();
			while ((textline = buffer.readLine()) != null) {
				//System.out.println(textline);
				if (textline.startsWith("%")) {
					if(name == null)
					{

					}
					else
						if (name.equalsIgnoreCase("%reference")) {
							sample.reference.sequences = sequences;
							sample.reference.sequenceNames = sequenceNames;
						} else if (name.equalsIgnoreCase("%mpd")) {
							sample.mpd.sequences = sequences;
							sample.mpd.sequenceNames = sequenceNames;
						}
						else
							if(name.equalsIgnoreCase("%posteriors"))
							{

							}
							else {
								Alignment alignment = new Alignment();
								alignment.sequences = sequences;
								alignment.sequenceNames = sequenceNames;
								sample.samples.add(alignment);
							}


					name = textline;
					if("%posteriors".equals(name))
					{
						for(int i = 0 ; (textline = buffer.readLine()) != null && !textline.startsWith("%") ; i++)
						{
							sample.posteriors.add(new Double(textline));
						}

						if(textline != null && textline.startsWith("%"))
						{
							name = textline;
							System.out.println("#"+name);

							//''sequences = new ArrayList<String>();
							//'sequenceNames = new ArrayList<String>();
						}
					}

					// some comment
					if(!alignmentString.equals(""))
					{
						parseAlignmentString(alignmentString, sequences, sequenceNames);
						alignmentString = "";

						sequences = new ArrayList<String>();
						sequenceNames = new ArrayList<String>();
					}
				} else {
					alignmentString += textline + "\n";
				}
			}


			if(name == null)
			{

			}
			else
				if (name.equalsIgnoreCase("%reference")) {
					sample.reference.sequences = sequences;
					sample.reference.sequenceNames = sequenceNames;
				} else if (name.equalsIgnoreCase("%mpd")) {
					sample.mpd.sequences = sequences;
					sample.mpd.sequenceNames = sequenceNames;
				}
				else
					if(name.equalsIgnoreCase("%posteriors"))
					{

					}
					else {
						Alignment alignment = new Alignment();
						alignment.sequences = sequences;
						alignment.sequenceNames = sequenceNames;
						sample.samples.add(alignment);
					}

			if(!alignmentString.equals(""))
			{
				parseAlignmentString(alignmentString, sequences, sequenceNames);
				alignmentString = "";

			}
			buffer.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		}

		return sample;
	}

	public static void parseAlignmentString(String alignmentString, ArrayList<String> sequences,  ArrayList<String> sequenceNames)
	{

		String seq = "";
		String name = "";
		String [] split = alignmentString.split("\n");
		for(int i = 0 ; i < split.length ; i++)
		{
			if(split[i].startsWith(">"))
			{
				name = split[i].substring(1);
				sequenceNames.add(name);
				if(!seq.equals(""))
				{
					sequences.add(seq);
					seq = "";
				}
				name = "";
			}
			else
			{
				seq += split[i].trim();
			}
		}
		if(!seq.equals(""))
		{
			sequences.add(seq);
		}
	}


	static class Alignment {
		ArrayList<String> sequences = new ArrayList<String>();
		ArrayList<String> sequenceNames = new ArrayList<String>();

		public String toString() {
			String ret = "";
			for (int i = 0; i < sequences.size(); i++) {
				ret += ">" + sequenceNames.get(i) + "\n";
				ret += sequences.get(i) + "\n";
			}
			return ret;
		}
	}	
}
