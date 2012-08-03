package statalign.distance;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import statalign.postprocess.utils.RNAFoldingTools;



public class Distance {
	final static String REFERENCE_FOLDER = "/home/ingolfur/Dropbox/RNA and StatAlign/TestRNAData/amap/ref";
	final static String MPD_FOLDER = "/home/ingolfur/Dropbox/RNA and StatAlign/TestRNAData/amap/mpd";
	final static String AMAP_FOLDER = "/home/ingolfur/Dropbox/RNA and StatAlign/TestRNAData/amap/amap";
	final static String AMAP4_FOLDER = "/home/ingolfur/Dropbox/RNA and StatAlign/TestRNAData/amap/amap4";
	final static String AMAP01_FOLDER = "/home/ingolfur/Dropbox/RNA and StatAlign/TestRNAData/amap/amap01";


	private static void ReadAlignments(ArrayList<Pair<ArrayList<String>,String>> alignments, String folder){
		File reffolder = new File(folder);
		File[] listOfRefFiles = reffolder.listFiles();
		for(File file  : listOfRefFiles){
			FileInputStream fstream;
			try {
				fstream = new FileInputStream(file.getAbsolutePath());
				DataInputStream in = new DataInputStream(fstream);
				BufferedReader br = new BufferedReader(new InputStreamReader(in));
				ArrayList<String> seq = new ArrayList<String>();
				String line = "";
				String nextline = "";
				br.readLine();
				while(nextline != null){
					nextline = br.readLine();
					line = "";
					while(nextline !=null && !nextline.contains(">") ){
						line += nextline;
						nextline = br.readLine();
					}
					seq.add(line.replace('.','-'));
				}

				String filename = file.getName().split(".fas")[0];
				filename = filename.replaceAll("yeah", "");
				filename = filename.replaceAll("removed", "");
				filename = filename.replaceAll("some4", "");
				filename = filename.replaceAll("some", "");
				alignments.add(new Pair<ArrayList<String>,String>(seq,filename));

				br.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}	

		}
	}



	public static void main(String[] args) {


		/*
		ArrayList<Pair<ArrayList<String>,String>> reference = new ArrayList<Pair<ArrayList<String>,String> > ();
		ReadAlignments(reference, REFERENCE_FOLDER);

		ArrayList<Pair<ArrayList<String>,String>> mpd = new ArrayList<Pair<ArrayList<String>,String> >();
		ReadAlignments(mpd, MPD_FOLDER);

		ArrayList<Pair<ArrayList<String>,String>> amap = new ArrayList<Pair<ArrayList<String>,String> >();
		ReadAlignments(amap, AMAP_FOLDER);

		ArrayList<Pair<ArrayList<String>,String>> amap4 = new ArrayList<Pair<ArrayList<String>,String> >();
		ReadAlignments(amap4, AMAP4_FOLDER);

		ArrayList<Pair<ArrayList<String>,String>> amap01 = new ArrayList<Pair<ArrayList<String>,String> >();
		ReadAlignments(amap01, AMAP01_FOLDER);

		Collections.sort(reference, new ComplexPair());
		Collections.sort(mpd, new ComplexPair());
		Collections.sort(amap, new ComplexPair());
		Collections.sort(amap4, new ComplexPair());
		Collections.sort(amap01, new ComplexPair());


		for(int i = 0; i<reference.size(); ++i){
			for(int j = 0; j<mpd.size(); ++j){
				for(int k = 0; k<amap.size(); ++k){
					for(int t = 0; t<amap01.size(); ++t){
						Pair<ArrayList<String>,String> refPair = reference.get(i);
						Pair<ArrayList<String>,String> mpdPair = mpd.get(j);
						Pair<ArrayList<String>,String> amapPair = amap.get(k);
						Pair<ArrayList<String>,String> amap01Pair = amap01.get(t);
						if(mpdPair.getRight().compareTo(refPair.getRight()) == 0 && amapPair.getRight().compareTo(refPair.getRight()) == 0 && amap01Pair.getRight().compareTo(refPair.getRight())==0){
							//System.out.println(Distance.AMA(refPair.getLeft(), mpdPair.getLeft()));
							//System.out.println(Distance.AMA(refPair.getLeft(), amapPair.getLeft()));
							//System.out.println(Distance.AMA(refPair.getLeft(), amap01Pair.getLeft()));
							//System.out.print("\t");
							System.out.println(amap01Pair.getRight());
							//System.out.println("----------------------------------");
						}
					}
				}
			}
		}
		 */
		ArrayList<String> seq = new ArrayList<String>();
		ArrayList<String> name = new ArrayList<String>();
		File folder = new File("/home/ingolfur/oxford_workspace/cmdStatAlign/seq5/seq5Clean/Cleaner/");
		File[] listOfFiles = folder.listFiles();
		Arrays.sort(listOfFiles);



		for(File i : listOfFiles){
			RNAFoldingTools.loadFastaSequences(i,seq , name);
			int length = 0;
			ArrayList<Integer> len = new ArrayList<Integer>();
			double lendist = 0;
			for(String s : seq){
				int seqlen = s.replace(".", "").length();
				len.add(seqlen);
				length += seqlen;
			}
			int count = 0;
			for(int j = 0; j < seq.size()-1; ++j ){
				for(int k = j+1; k < seq.size(); ++k ){
					lendist += Math.abs(len.get(j) -len.get(k)) / (double)(Math.max(len.get(j), len.get(k)) -1);
					count++;
				}
			}
			lendist = lendist / (double) count;




			System.out.println(i.getName());
			//System.out.println(length);
			//System.out.println(lendist);
			//System.out.println(Distance.sequenceSimilarityScore(seq));
			//System.out.println("-----");
			seq.clear();
		}

		//System.out.println(sequenceSimilarityScore(referenceA1));





		/*
		ArrayList<String> referenceA1 = new ArrayList<String>();
		ArrayList<String> mpd1 = new ArrayList<String>();
		ArrayList<String> amapnormTest1 = new ArrayList<String>();
		ArrayList<String> amapseq4Test1 = new ArrayList<String>();

		ArrayList<String> mpd17 = new ArrayList<String>();
		ArrayList<String> referenceA17 = new ArrayList<String>();
		ArrayList<String> amapnormTest17 = new ArrayList<String>();
		ArrayList<String> amapseq4Test17 = new ArrayList<String>();

		ArrayList<String> mpd14 = new ArrayList<String>();
		ArrayList<String> referenceA14 = new ArrayList<String>();
		ArrayList<String> amapnormTest14 = new ArrayList<String>();
		ArrayList<String> amapseq4Test14 = new ArrayList<String>();

		ArrayList<Pair<Double,Integer>> refAndMPD = new ArrayList<Pair<Double,Integer>>();
		ArrayList<Pair<Double,Integer>> refAndAMAP = new ArrayList<Pair<Double,Integer>>();
		ArrayList<Pair<Double,Integer>> refAndAMAP4 = new ArrayList<Pair<Double,Integer>>();



		mpd1.add("GGGUGCUUGAGGCUGUCUGCCUCGGGCA--UGCCAC---UGUAAGGCAGACAGAGAAAAGCCCCAGUU-AACAUUACGCGUCCUGCAAGACGCUUAACAUUAA-UCUGAGGCCCAAUCU-AUGUCUCA-CAAA---UGU---AGGUUAGCCUCUUACGUGCCGAAAGGCAA----GGAGAAGCAGGCU-AUGAAG");
		mpd1.add("GGGUGCUUGAAGCUGUCUGCUUUAAGUGCUUGCAUCAGGCUGAGAGUAGGCAGAGAAAAGCCCCGUAUCAAUGUUA----------------------AUCAA-UACGAGGCCCUCUGUAAUGCACGA-CAAC---AUU---ACGGUAGCCUUUUACCCGCCGAAAGGCAA----GGAG------GCUGAAGAUG");
		mpd1.add("GGGUGCUUGAGACUGUUUGUCUCAGGUAUU---UAC---CAAAAGGCAGACAGAGAAAAGCCCCACCUGACUAUAA----------------------AUCAAAAGUGCAUUGCAC-CC-AUUAUUGAUCUCUUCAAUAACGAAGCUAUCCCCU-----ACAGUAUUUCAG----AACG------UCCAACCAUG");
		mpd1.add("GGGUGCUUGAGGCUGUCUGCCUCGGGCA--UGCCAC---CGUAAGGCAGACAGAGAAAAGCCCCAGUU-AACAUUACGCGUCCUGCAAGACGCCUAACAUUAA-UCUGAGGCCAAUUUC-AUGCUAGA-CACA---UGU---AGGUUAGCCUCUUACGCGCCGAAAGGCAA----GGAGAAGCA-GCU-AUG---");
		mpd1.add("GGGUGCUUGAGACUGUUUGUCUCAGGUAUU---CAC---CGAAAGGCAGACAGAGAAAAGCCCCACCUGACUAUAA----------------------AUCAA-AGUGAGGCU-ACCCU-AUGCCUGAACACC---AUA---AGGUUAGCCUCUUACUCGUUGGAAAUCAACACAGGGG------GCUGGGAAUG");
		refAndMPD.add(new Pair<Double,Integer>(Distance.AMA(mpd1, referenceA1),1));

		amapnormTest1.add("GGGUGCUUGAAGCUG----UCUGCUUUAAGUGCUUGCAUCAGGC--------------UGAGAGU-----AGGCAGAGAAAAGCCCCGUAUCAAUGUUAAUCAAU----ACG------------------------AGG--------CCCUCUGUAAUG--CACGAC----AACAUU----------------------------------------------------A--CGG-UAGCCU----UUUAC----------------------------------CCGCCGAAAGGC------------AAGGA------GGCU---GAAGAUG");
		amapnormTest1.add("GGGUGCUUGAGACUGUUUGUCU-----------------CAGG-UAUUUA--CCAA---------AAGGCAGACAGAGAAAAGCCCCA---------------------------CCUGACUAUAAAUCAAAAGUG---CAUUGC----------------------AC------------------CCAUUAUUG-A--U--CUCUUCA-----------------------------------------AUAACGAAGCUAUCCCCUACAGUAUUUCA-----------------GAACGUCCAACC----------------------AUG");
		amapnormTest1.add("GGGUGCUUGAGGCUG----UCUG-------------CCUCGGGC-----AUGCCACUGU------AAGGCAGACAGAGAAAAGCCCCA-------GUUAA-----CAUUACGCGUCC----------------------------------------UG--CAAGACGCUUAACAUUAAUCUGAGGCCC-------AA--U------------CUAUGUCUCACAA-AU-GUAGGUUAGCCUCUUA----CGU----------------------------------GCCGAAAGGC------------AAGGAGAAGCAGGC-UAUGAAG---");
		amapnormTest1.add("GGGUGCUUGAGACUGUUUGUCU-----------------CAGG-UAUUCA--CCGA---------AAGGCAGACAGAGAAAAGCCCCA---------------------------CCUGACUAUAAAUCAA-AGUGAGG-----CUACCCUAUG-CC--UGAAC-------A-----------------------------------------------------CCAUA--AGGUUAGCCUCUUA----C-----------------------------UCGUU------------GGAAAUCAACACAGG-G------GGCU---GGGAAUG");
		amapnormTest1.add("GGGUGCUUGAGGCUG----UCUG-------------CCUCGGGC-----AUGCCACCGU------AAGGCAGACAGAGAAAAGCCCCA-------GUUAA-----CAUUACGCGUCC----------------------------------------UG--CAAGACGCCUAACAUUAAUCUGAGGCCA-------A-UU-UC-------AUGCUA-----GACAC-AU-GUAGGUUAGCCUCUUA----C----------------------------------GCGCCGAAAGGC------------AAGGAGAAGCA-GC-UAU------G");
		refAndAMAP.add(new Pair<Double,Integer>(Distance.AMA(amapnormTest1, referenceA1),1));

		amapseq4Test1.add("GGGUGCUUGAAGCUG-----UCUGCUUUAAGUGCUUGC--AU-CAGGCUGAGA----------GU--------------AGGCAGAGAAAAGCCCCGUAUCAAU-GUUAAUCAAUA---------CGA--------------------------GGC-----CCUCUGUAAUGCACGAC----------------------------------------------------------------------------AAC----AUUACGGU--------------------------------------------------AGCC------UUUUAC-----C--CGCCGAAAGGCAAGG-----------------AGGCU--------------GAAGAU-G");
		amapseq4Test1.add("GGGUGCUUGAGACU-GUUUGUC--------------------UC---------AGGUAUUUA-----CCAA---AAGGCAGACAGAGAAAAGCCCC--------A--------------------------CCUGACUAUAAAUCAA-AAG-----------------------------------UGCAUUGCACCCAUUAUUGAUCUCUUC-AAUAACGAAGCU-AUCCCCUACAGUAUUUCAGAACGU------------------------------------------------------------------------------------------------------------------CCAACC----------------------------AU-G");
		amapseq4Test1.add("GGGUGCUUGAGGCUG-----UCU---------------GCCU-CGGGC--------------A--UGCCACUGUAAGGCAGACAGAGAAAAGCCCC--------AG----------UUAACAUUACG-CGUCC--------------------------------------UGCAAGACGCUUA------------------------------------------------------------------------AC----AUU-----AAUCUGAGGCCCAAUCUAUGUCUC-------------ACAAAUGUAGGUUAGCCUCUUAC------------GU-GCCGAAAGGCAAGG----------------------AGAAGCAGG-CUAUGAA---G-");
		amapseq4Test1.add("GGGUGCUUGAGACU-GUUUGUC--------------------UC---------AGGUAUUCA-----CCGA---AAGGCAGACAGAGAAAAGCCCC--------A--------------------------CCUGACUAUAAAUCAAA--GUGAG--GCUACCCUAU-----------------GC---------------------------C------------U------------------------GAACA--CCAU---AA-------------------------------------------------GGUUAGCCUCUUAC------UCGUU------------------GGAAAUCAAC--ACAGGGGGCU--------------GGGAAU-G");
		amapseq4Test1.add("GGGUGCUUGAGGCUG-----UCU---------------GCCU-CGGGC--------------A--UGCCACCGUAAGGCAGACAGAGAAAAGCCCC--------AG----------UUAACAUUACG-CGUCC--------------------------------------UGCAAGACGCCUA------------------------------------------------------------------------AC----AUU-----AAUCUGAGGCC-------------AAUUUCAUGCUAGACACAUGUAGGUUAGCCUCUUAC-----------G--CGCCGAAAGGCAAGG----------------------AGAAGCA--GCUAU-------G");
		refAndAMAP4.add(new Pair<Double,Integer>(Distance.AMA(amapseq4Test1, referenceA1),1));

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		referenceA17.add("C--------------------------------------------G-U-G-----------C--------------------------------------------AGCCUCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCUGGGGCGACCGGGUCCUUUCUUGGAACAACCCGCUCAAUACCCAGAGAUUUGGGCGUGCCCCCGCAAGACCACUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCCUAAACCUCA");
		referenceA17.add("---------------------------------------------G-GCG-----------UUA---------------------------GUAUGAGUGUUGUACAGCCUCCAGGCCCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCCGGGAUGACCGGGUCCUUUCUUGGAUUAACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCAAGACUGCUAG-------------------------CCGAGUAGUGUU--------------GG----GU--CGCGAAAGG-----------------------------------------");
		referenceA17.add("GCCAGCCCCUAAUGGGGCGACACUCCACCAUGAUCACUCCCCUGUGAG-GAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUGAGUGUCGUGCAGCCUCCAGGACUCCCCCUCCCGGGAGAGCCAUAGUAGUCUGCGGAACCGGUGAGUACACCGGAAUUGCCAGGAUGACCGGGUCCUUUCUUGGAUCAACCCGCUCGAUGCCUGGAGAUUUGGGCGUGCCCCCGCAAGAUUGCUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCCAAAACCCCA");
		referenceA17.add("---------------------------------------------G-GCG-----------U---------------------U------AGUAUAUAGUGCGUGCAGCCUCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUUGCCAGGAUGACCGGGUCCUUUCUUGGAUUAACCCGCUCAGUGCCUGGAGAUUUGGGCGUGCCCCCGCGAGACUGCUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAUCAUGAGCACACUUCCAAAACCCCA");
		referenceA17.add("---------------------------------------------------------------------------------------------------------CAGCCCCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCCGGGAUGACCGGGUCCUUUCUUGGAUAAACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCGAGACCGUUAGCCGAGUAGUGUUGGGUCGCAAAAGGC---------CU---------------------U-------GUGG-----------------------------------------");

		mpd17.add("---------------------------------------------------------------------------------------------------------CAGCCCCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCCGGGAUGACCGGGUCCUUUCUUGGAUAAACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCGAGACCGUUAGCCGAGUAGUGUUGGGUCGCAAAAGGCCUUGUGG------------------------------------------------------------------------------");
		mpd17.add("GCCAGCCCCUAAUGGGGCGACACUCCACCAUGAUCACUCCCCUGUGAGGAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUG-AGUGUCGUGCAGCCUCCAGGACUCCCCCUCCCGGGAGAGCCAUAGUAGUCUGCGGAACCGGUGAGUACACCGGAAUUGCCAGGAUGACCGGGUCCUUUCUUGGAUCAACCCGCUCGAUGCCUGGAGAUUUGGGCGUGCCCCCGCAAGAUUGCUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCCAAAACCCCA");
		mpd17.add("-----------------------------------------------------------------------------------GGCGUUAGUAUAUAGUG-CGUGCAGCCUCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUUGCCAGGAUGACCGGGUCCUUUCUUGGAUUAACCCGCUCAGUGCCUGGAGAUUUGGGCGUGCCCCCGCGAGACUGCUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAUCAUGAGCACACUUCCAAAACCCCA");
		mpd17.add("-----------------------------------------------------------------------------------------------------CGUGCAGCCUCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCUGGGGCGACCGGGUCCUUUCUUGGAACAACCCGCUCAAUACCCAGAGAUUUGGGCGUGCCCCCGCAAGACCACUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCCUAAACCUCA");
		mpd17.add("-----------------------------------------------------------------------------------GGCGUUAGUAUG-AGUGUUGUACAGCCUCCAGGCCCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCCGGGAUGACCGGGUCCUUUCUUGGAUUAACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCAAGACUGCUAGCCGAGUAGUGUUGGGUCGCGAAAGG--------------------------------------------------------------------------------------");
		refAndMPD.add(new Pair<Double,Integer>(Distance.AMA(mpd17, referenceA17),17));

		amapnormTest17.add("------------------------------------------------------------------------------------------------------CGUGCAGCCUCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCUGGGGC-GACCGGGUCCUUUCUUGGAAC--AACCCGCUCAAUACCCAGAGAUUUGGGCGUGCCCCCGCAAGAC--CACUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCCUAAA--CCUCA");
		amapnormTest17.add("-----------------------------------------------------------------------------------GGCGUUAGUAUG-AGUGUU-GUACAGCCUCCAGGCCCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGC-CGGGAUGACCGGGUCCUUUCUUGGA--UUAACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCAAGACUGC--UAGCCGAGUAGUGUUGGGUCGCGAAAGG----------------------------------------------------------------------------------------");
		amapnormTest17.add("GCCAGCCCCUAAUGGGGCGACACUCCACCAUGAUCACUCCCCUGUGAGGAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUG-AGUGU-CGUGCAGCCUCCAGGACUCCCCCUCCCGGGAGAGCCAUAGUAGUCUGCGGAACCGGUGAGUACACCGGAAUUGC-CAGGAUGACCGGGUCCUUUCUUGGA--UCAACCCGCUCGAUGCCUGGAGAUUUGGGCGUGCCCCCGCAAGAUUGC--UAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCC-AAAACCC-CA");
		amapnormTest17.add("-----------------------------------------------------------------------------------GGCGUUAGUAUAUAGUG--CGUGCAGCCUCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUUGC-CAGGAUGACCGGGUCCUUUCUUGGA--UUAACCCGCUCAGUGCCUGGAGAUUUGGGCGUGCCCCCGCGAGACUGC--UAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAUCAUGAGCACACUUCC-AAAACCC-CA");
		amapnormTest17.add("--CAGCCC--------------------------------------------------------------------------------------------------------CCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGC-CGGGAUGACCGGGUCCUUUCUUGGA--UAAACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCGAGAC--CGUUAGCCGAGUAGUGUUGGGUCGCAAAAGGCCUUGUGG--------------------------------------------------------------------------------");
		refAndAMAP.add(new Pair<Double,Integer>(Distance.AMA(amapnormTest17, referenceA17),17));

		amapseq4Test17.add("-------------------------------------------------------------------------------------------------------------CGUGCAGCCUCCAGGAC-CCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCUGG-----GGC-GACCGGGUCCUUUCUUGGAAC--A-ACCCGCUCAAUACCCAGAGAUUUGGGCGUGCCCCCGCAAGACCAC----UAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCCUAAA---CCUCA");
		amapseq4Test17.add("-----------------------------------------------------------------------------------GGCGUUAGUAU-----GAGU-GUUGU---ACAGCCUCCAGGCC-CCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCG------CCGGGAUGACCGGGUCCUUUCUUGGA--UU-AACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCAAGAC----UGCUAGCCGAGUAGUGUUGGGUCGCGAAAGG-----------------------------------------------------------------------------------------");
		amapseq4Test17.add("GCCAGCCCCUAAUGGGGCGACACUCCACCAUGAUCACUCCCCUGUGAGGAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAU-----GAGUG-U---CGUGCAGCCUCCAGGACU-CCCCCUCCCGGGAGAGCCAUAGUAGUCUGCGGAACCGGUGAGUACACCGGAAU------UGCCAGGAUGACCGGGUCCUUUCUUGGA--UCA-ACCCGCUCGAUGCCUGGAGAUUUGGGCGUGCCCCCGCAAGA----UUGCUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCC-AA-AACCC-CA");
		amapseq4Test17.add("-----------------------------------------------------------------------------------GGCGUUAGUAUAUAGU----G-----CGUGCAGCCUCCAGGAC-CCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAU------UGCCAGGAUGACCGGGUCCUUUCUUGGA--UU-AACCCGCUCAGUGCCUGGAGAUUUGGGCGUGCCCCCGCGAGAC----UGCUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAUCAUGAGCACACUUCC-AA-AACCC-CA");
		amapseq4Test17.add("--CAGCCC---------------------------------------------------------------------------------------------------------------CCAGGAC-CCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCG------CCGGGAUGACCGGGUCCUUUCUUGGA--UAA-ACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCGAGACCGU----UAGCCGAGUAGUGUUGGGUCGCAAAAGGCCUUGUGG---------------------------------------------------------------------------------");
		refAndAMAP4.add(new Pair<Double,Integer>(Distance.AMA(amapseq4Test17, referenceA17),17));

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		referenceA14.add("GCUUUCUCCCAGGUG-GGACUGGGGGUAAAGCGAGCGUCGAACUGGUAGCCGCCGGUUCGGCGCUCAU--UUUUUCGU-UU");
		referenceA14.add("GCUUUCUCCCAGGUGUGUGCUGGGUGAUAAGCGAAAGUC-AUCGGGUUGCCGCCCGGU-GGCUUUCUUCGUUUUUCAUUGU");
		referenceA14.add("GCUUUCUCCCAGCUGAUGACUGGGGG-UUAGCCGACGCCUCGACAGUUACAGCUGUCGGGGCGUUCAACAUUUUUCAA-GU");
		referenceA14.add("GCUUUCUCCCAGCUGAUGACUGGGGG-UUAGCCGACGCCCUGUGAGUUCCCGCUCACGGGGCGUUCAAC-UUUUUCAG-GU");
		referenceA14.add("GCUUUCUCCCAGGUAGUGGCUGGGUGUAAAGCGAGCGUCGAACUGGUUGCAGCCGGUUCGGCGCUCAUU-UUUUUCGU-UU");

		mpd14.add("GCUUUCUCCCAGCUGAUGACUGGGGGUUAGCCGA-CGCC-UCGACAGUUACAGCUGUCGGGGCGUUCAACAUUUUUCAA-GU");
		mpd14.add("GCUUUCUCCCAGGUGG-GACUGGGGGUAAAGCGAGCGUCGAACUG-GUAGCCGCCGGUUCGGCGCUCAUU--UUUUCGU-UU");
		mpd14.add("GCUUUCUCCCAGGUGUGUGCUGGGUGAUAAGCGAAAGUC-AUCGG-GUUGCCGCCCGGU-GGCUUUCUUCGUUUUUCAUUGU");
		mpd14.add("CUUUCUCCCAGGUAGUGGCUGGGUGUAAAGCGAGCGUCGAACUG-GUUGCAGCCGGUUCGGCGCUCAUU-UUUUUCGU-UU");
		mpd14.add("GCUUUCUCCCAGCUGAUGACUGGGGGUUAGCCGA-CGCC-CUGUGAGUUCCCGCUCACGGGGCGUUCAAC-UUUUUCAG-GU");
		refAndMPD.add(new Pair<Double,Integer>(Distance.AMA(mpd14, referenceA14),14));

		amapnormTest14.add("GCUUUCUCCCAGGUGG-G--ACUGGG--GG-UAA--AGCGA-GC---GUCGAACUG----GUA-----GCCGC---C-GGUUCGG----CG-------------------CUCAUUU-UUUCGUUU--------");
		amapnormTest14.add("GCUUUCUCCCAGGUG-U-GUGCUGGGUG-A-U---AAGCGA-AAGUCAUCG-------------GGUUGCCGC---CCG-----GUGGC---------------UUUC-----------UUCGUUUUUCAUUGU");
		amapnormTest14.add("GCUUUCUCCCAGCUGAUG--ACUGGG--GGU---UAGCCGACG---CCUCGA----C-------AGUUACAGCUGUC-GG------GGC-GUUCAACAUUUUUC----AA----------------------GU");
		amapnormTest14.add("GCUUUCUCCCAGCUGAUG--ACUGGG--GGU---UAGCCGACG---C-----C---CUGU---GAGUUCCCGCUCA-CGG------GGC-GUUCAAC-UUUUUC----AG----------------------GU");
		amapnormTest14.add("GCUUUCUCCCAGGUAGUG--GCUGGGUG---UAA--AGCGA-GC---GUCGAACUG---------GUUGCAGC---C-GGUUCGG----CG-------------------CUCAUUUUUUUCGUUU--------");
		refAndAMAP.add(new Pair<Double,Integer>(Distance.AMA(amapnormTest14, referenceA14),14));

		amapseq4Test14.add("GCUUUCUCCCAGGUGG-GA-----CUGGG-----GGUAAAGCGAGC--------------GUCGAACU-G-----GUA------GCCGC----C-------------------------------------GGUUCGGCGCUCAUUU-UUUCGUU---------U");
		amapseq4Test14.add("GCUUUCUCCCAGGU-----GUGUGCUGGGUGAUA-----AGCGAAAGUC-----------AUCG----G-----------GUUG-CCGC----C--CGGUG----------------------GCUUUCUU--------------------CGUUUUUCAUUGU-");
		amapseq4Test14.add("GCUUUCUCCCAGCUGAUGA-----CUGGG-----GG-------------UUAGCCGACGCCUCGA------C-------AGUUA-CAGCUGUC--------GGGGCGUUCAACAUUUUUCAAG----------------------------------------U-");
		amapseq4Test14.add("GCUUUCUCCCAGCUGAUGA-----CUGGG-----GG-------------UUAGCCGACGC----------CCUGU---GAGUUC-CCGCU----CAC----GGGGCGUUCAAC-UUUUUCAGG----------------------------------------U-");
		amapseq4Test14.add("GCUUUCUCCCAGGUAGUGG-----CUGGGUG-----UAAAGCGAGC--------------GUCGAACU-G----------GUUG-CAGC----C-------------------------------------GGUUCGGCGCUCAUUUUUUUCGUU---------U");
		refAndAMAP4.add(new Pair<Double,Integer>(Distance.AMA(amapseq4Test14, referenceA14),14));





		for(int i = 0; i<refAndMPD.size(); ++i){
			System.out.println("Test" + refAndMPD.get(i).getRight());
			System.out.println("The AMA between the reference alignment and the MPD from Test" +refAndMPD.get(i).getRight() + ": " + refAndMPD.get(i).getLeft());
			System.out.println("The AMA between the reference alignment and the AMAP from Test" +refAndAMAP.get(i).getRight() + ": " + refAndAMAP.get(i).getLeft());
			System.out.println("The AMA between the reference alignment and the AMAP4 from Test" +refAndAMAP4.get(i).getRight() + ": " + refAndAMAP4.get(i).getLeft());
			System.out.println();
		}

		System.out.println("The distance between the reference alignment and the MPD from Test1: "+Distance.multiDistance(referenceA1,mpd1));
		//System.out.println(Distance.distance(A1.toArray(),A2.toArray()));
		System.out.println("The AMA between the reference alignment and the MPD from Test1: "+Distance.AMA(referenceA1,mpd1));
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println("-----------------------------------------------------------------------------------\n");

		System.out.println("The distance between the reference alignment and the AMAP-4 from Test1: "+Distance.multiDistance(referenceA1,amapseq4Test1));
		//System.out.println(Distance.distance(A1.toArray(),A2.toArray()));
		System.out.println("AMA between those alignments: "+Distance.AMA(referenceA1,amapseq4Test1));
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println("The distance between the Normal AMAP and the reference alignment from Test1: "+Distance.multiDistance(amapnormTest1,referenceA1));
		//System.out.println(Distance.distance(A1.toArray(),A2.toArray()));
		System.out.println("AMA between those alignments: "+Distance.AMA(amapnormTest1,referenceA1));
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println();
		System.out.println();

		System.out.println("The distance between the reference alignment and the MPD from Test17: "+Distance.multiDistance(referenceA17,mpd17));
		//System.out.println(Distance.distance(A1.toArray(),A2.toArray()));
		System.out.println("AMA between those alignments: "+Distance.AMA(referenceA17,mpd17));
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println("-----------------------------------------------------------------------------------\n");
		System.out.println("The distance between the Normal AMAP and the reference alignment from Test17: "+Distance.multiDistance(amapnormTest17,referenceA17));
		System.out.println("AMA between those alignments: "+Distance.AMA(amapnormTest17,referenceA17));
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println("The distance between the  AMAP-4 and the reference alignment from Test17: "+Distance.multiDistance(amapseq4Test17,referenceA17));
		System.out.println("AMA between those alignments: "+Distance.AMA(amapseq4Test17,referenceA17));
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println();
		System.out.println();


		System.out.println("The distance between the reference alignment and the MPD from Test17: "+Distance.multiDistance(referenceA14,mpd14));
		//System.out.println(Distance.distance(A1.toArray(),A2.toArray()));
		System.out.println("AMA between those alignments: "+Distance.AMA(referenceA14,mpd14));
		System.out.println("-----------------------------------------------------------------------------------");
		System.out.println("-----------------------------------------------------------------------------------\n");
		System.out.println("The distance between the Normal AMAP and the reference alignment from Test17: "+Distance.multiDistance(amapnormTest14,referenceA14));
		System.out.println("AMA between those alignments: "+Distance.AMA(amapnormTest14,referenceA14));
		System.out.println("-----------------------------------------------------------------------------------");




		ArrayList<String> B1 = new ArrayList<String>();
		ArrayList<String> B2 = new ArrayList<String>();

		B1.add("C--------------------------------------------G-U-G-----------C--------------------------------------------AGCCUCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCUGGGGCGACCGGGUCCUUUCUUGGAACAACCCGCUCAAUACCCAGAGAUUUGGGCGUGCCCCCGCAAGACCACUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCCUAAACCUCA");
		B1.add("---------------------------------------------G-GCG-----------UUA---------------------------GUAUGAGUGUUGUACAGCCUCCAGGCCCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCCGGGAUGACCGGGUCCUUUCUUGGAUUAACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCAAGACUGCUAG-------------------------CCGAGUAGUGUU--------------GG----GU--CGCGAAAGG-----------------------------------------");
		B1.add("GCCAGCCCCUAAUGGGGCGACACUCCACCAUGAUCACUCCCCUGUGAG-GAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUGAGUGUCGUGCAGCCUCCAGGACUCCCCCUCCCGGGAGAGCCAUAGUAGUCUGCGGAACCGGUGAGUACACCGGAAUUGCCAGGAUGACCGGGUCCUUUCUUGGAUCAACCCGCUCGAUGCCUGGAGAUUUGGGCGUGCCCCCGCAAGAUUGCUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCCAAAACCCCA");	
		B1.add("---------------------------------------------G-GCG-----------U---------------------U------AGUAUAUAGUGCGUGCAGCCUCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUUGCCAGGAUGACCGGGUCCUUUCUUGGAUUAACCCGCUCAGUGCCUGGAGAUUUGGGCGUGCCCCCGCGAGACUGCUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAUCAUGAGCACACUUCCAAAACCCCA");
		B1.add("---------------------------------------------------------------------------------------------------------CAGCCCCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCCGGGAUGACCGGGUCCUUUCUUGGAUAAACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCGAGACCGUUAGCCGAGUAGUGUUGGGUCGCAAAAGGC---------CU---------------------U-------GUGG-----------------------------------------");

		B2.add("-----------------------------------------------------------------------------------------------------CGUGCAGCCUCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCUGGGGCGACCGGGUCCUUUCUUGGAACAACCCGCUCAAUACCCAGAGAUUUGGGCGUGCCCCCGCAAGACCACUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCCUAAACCUCA");
		B2.add("-----------------------------------------------------------------------------------GGCGUUAGUAUAUAGUG-CGUGCAGCCUCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUUGCCAGGAUGACCGGGUCCUUUCUUGGAUUAACCCGCUCAGUGCCUGGAGAUUUGGGCGUGCCCCCGCGAGACUGCUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAUCAUGAGCACACUUCCAAAACCCCA");
		B2.add("GCCAGCCCCUAAUGGGGCGACACUCCACCAUGAUCACUCCCCUGUGAGGAACUACUGUCUUCACGCAGAAAGCGUCUAGCCAUGGCGUUAGUAUG-AGUGUCGUGCAGCCUCCAGGACUCCCCCUCCCGGGAGAGCCAUAGUAGUCUGCGGAACCGGUGAGUACACCGGAAUUGCCAGGAUGACCGGGUCCUUUCUUGGAUCAACCCGCUCGAUGCCUGGAGAUUUGGGCGUGCCCCCGCAAGAUUGCUAGCCGAGUAGUGUUGGGUCGCGAAAGGCCUUGUGGUACUGCCUGAUAGGGUGCUUGCGAGUGCCCCGGGAGGUCUCGUAGACCGUGCAACAUGAGCACACUUCCAAAACCCCA");
		B2.add("-----------------------------------------------------------------------------------GGCGUUAGUAUG-AGUGUUGUACAGCCUCCAGGCCCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCCGGGAUGACCGGGUCCUUUCUUGGAUUAACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCAAGACUGCUAGCCGAGUAGUGUUGGGUCGCGAAAGG--------------------------------------------------------------------------------------");
		B2.add("---------------------------------------------------------------------------------------------------------CAGCCCCCAGGACCCCCCCUCCCGGGAGAGCCAUAGUGGUCUGCGGAACCGGUGAGUACACCGGAAUCGCCGGGAUGACCGGGUCCUUUCUUGGAUAAACCCGCUCAAUGCCCGGAAAUUUGGGCGUGCCCCCGCGAGACCGUUAGCCGAGUAGUGUUGGGUCGCAAAAGGCCUUGUGG------------------------------------------------------------------------------");
		 */
	}


	/**
	 * Given a set of sequences and an AMA score, returns the corresponding distance.
	 * @param sequences
	 * @param amaScore
	 * @return
	 * @param sim
	 * @return 
	 */
	public static double amaScoreToMultiDistance(ArrayList<String> sequences, double amaScore)
	{
		double lengthSum = 0;
		for(int i = 0 ; i < sequences.size() ; i++)
		{
			lengthSum += sequences.get(i).replaceAll("-", "").length();
		}
		double ksub1 = sequences.size() - 1;

		return (1-amaScore)*ksub1*lengthSum;
	}

	public static int multiDistance(String [] A, String [] B){
		ArrayList<String> arListA = new ArrayList <String> (Arrays.asList(A));
		ArrayList<String> arListB = new ArrayList <String> (Arrays.asList(B));
		return Distance.multiDistance(arListA, arListB);

	}

	public static int multiDistance(ArrayList<String>  A, ArrayList<String>  B){
		ArrayList<String> cloneA = new ArrayList<String>();
		ArrayList<String> cloneB = new ArrayList<String>();

		for(String i : A){
			cloneA.add(i);
		}

		for(String i : B){
			cloneB.add(i);
		}

		sortSeq(cloneA,cloneB);
		String [] tempA = new String[2]; 
		String [] tempB = new String[2];
		int d = 0;
		int k = A.size(); //This says how many sequences the alignments have
		for(int i=0; i<k-1; ++i){
			for(int j=i+1; j<k; ++j){
				tempA[0] = A.get(i);
				tempA[1] = A.get(j);
				tempB[0] = B.get(i);
				tempB[1] = B.get(j);
				d += Distance.distance(tempA, tempB);
			}
		}

		return d;
	}

	private static void sortSeq(ArrayList<String>  A, ArrayList<String>  B){

		List< Pair < String  , Integer> > seq1 = new ArrayList< Pair < String , Integer> >();
		List< Pair < String  , Integer> > seq2 = new ArrayList< Pair < String  , Integer> >();
		ArrayList<String> newA = new ArrayList<String>();
		ArrayList<String> newB = new ArrayList<String>();
		for(int i = 0; i<A.size(); ++i){
			seq1.add( new Pair<String,Integer >  (A.get(i).replaceAll("-",""),i));
			seq2.add( new Pair<String,Integer >  (B.get(i).replaceAll("-",""),i));
		}
		Collections.sort(seq1, new PairCompare());
		Collections.sort(seq2, new PairCompare());

		for(int i = 0; i<A.size(); ++i){
			int index1 =seq1.get(i).getRight();
			int index2 =seq2.get(i).getRight();

			newA.add(A.get(index1));
			newB.add(B.get(index2));
		}
		for(int i = 0; i<A.size(); ++i){
			A.set(i, newA.get(i));
			B.set(i, newB.get(i));
		}
	}

	private static int distance(String[] A, String[] B){
		String seqA1 = A[0].replaceAll("-", "");
		String seqA2 = A[1].replaceAll("-", "");

		int homoSize = Distance.HomoUnionSize(A, B);
		int delandInsertSize = deletionUnionAndInsertionUnionSize(A,B); 

		int dist = seqA1.length() + seqA2.length() - 2 * homoSize - delandInsertSize;
		return dist;

	}

	private static int HomoUnionSize(String[] A, String[] B) {
		ArrayList<Pair<Integer,Integer>> HomoForA = new ArrayList<Pair<Integer,Integer>>();
		int seqIndex1 = 0;
		int seqIndex2 = 0;
		char seq1Char;
		char seq2Char;
		for(int i =0; i<A[0].length(); ++i ){
			seq1Char = A[0].charAt(i);
			seq2Char = A[1].charAt(i);

			if(seq2Char != '-' && seq1Char != '-') {
				HomoForA.add( new Pair<Integer,Integer>(seqIndex1,seqIndex2));
			}

			if(seq1Char != '-'){
				seqIndex1++;
			}


			if(seq2Char != '-'){
				seqIndex2++;
			}

		}

		ArrayList<Pair<Integer,Integer>> HomoForB = new ArrayList<Pair<Integer,Integer>>();
		seqIndex1 = 0;
		seqIndex2 = 0;
		for(int i =0; i<B[0].length(); ++i ){
			seq1Char = B[0].charAt(i);
			seq2Char = B[1].charAt(i);

			if(seq2Char != '-' && seq1Char != '-') {
				HomoForB.add( new Pair<Integer,Integer>(seqIndex1,seqIndex2));
			}

			if(seq1Char != '-'){
				seqIndex1++;
			}

			if(seq2Char != '-'){
				seqIndex2++;
			}
		}

		ArrayList<Pair<Integer,Integer>> HomoUnion = new ArrayList<Pair<Integer,Integer>>();
		for(int j=0;j<HomoForA.size();++j)
		{	
			Pair<Integer,Integer> common = HomoForA.get(j);
			if(HomoForB.contains(common)){
				HomoUnion.add(common);
			}
		}
		return HomoUnion.size();

	}

	private static int deletionUnionAndInsertionUnionSize(String[] A, String[] B) {
		ArrayList<Integer> deletionForA = new ArrayList<Integer>();
		ArrayList<Integer> insertionForA = new ArrayList<Integer>();
		char seq1Char;
		char seq2Char;
		int seqIndex1 = 0;
		int seqIndex2 = 0;
		for(int i =0; i<A[0].length(); ++i ){
			seq1Char = A[0].charAt(i);
			seq2Char = A[1].charAt(i);

			if(seq1Char == '-' && seq2Char == '-'){
				continue;
			}

			if(seq1Char == '-'){
				insertionForA.add(seqIndex2);
			}
			else{
				seqIndex1++;
			}

			if(seq2Char == '-'){
				deletionForA.add(seqIndex1-1);
			}
			else{
				seqIndex2++;
			}


		}

		seqIndex1=0;
		seqIndex2=0;
		ArrayList<Integer> deletionForB = new ArrayList<Integer>();
		ArrayList<Integer> insertionForB = new ArrayList<Integer>();
		for(int i =0; i<B[0].length(); ++i ){
			seq1Char = B[0].charAt(i);
			seq2Char = B[1].charAt(i);

			if(seq1Char == '-' && seq2Char == '-'){
				continue;
			}


			if(seq2Char == '-'){
				deletionForB.add(seqIndex1);
			}
			else{
				seqIndex2++;
			}

			if(seq1Char == '-'){
				insertionForB.add(seqIndex2-1);
			}
			else{
				seqIndex1++;
			}
		}



		ArrayList<Integer> deletionUnion = new ArrayList<Integer>();
		for(int j=0;j<deletionForA.size();++j)
		{	
			int common = deletionForA.get(j);
			if(deletionForB.contains(common)){
				deletionUnion.add(common);
			}
		}

		ArrayList<Integer> insertionUnion = new ArrayList<Integer>();
		for(int j=0;j<insertionForA.size();++j)
		{	
			int common = insertionForA.get(j);
			if(insertionForB.contains(common)){
				insertionUnion.add(common);
			}
		}
		return deletionUnion.size() + insertionUnion.size();
	}

	public static double AMA(String[] A, String[] B){
		ArrayList<String> arListA = new ArrayList <String> (Arrays.asList(A));
		ArrayList<String> arListB = new ArrayList <String> (Arrays.asList(B));
		return Distance.AMA(arListA, arListB);
	}

	public static double AMA(ArrayList<String>  A, ArrayList<String>  B){
		int sum = 0;
		for(int i = 0; i<A.size(); ++i){
			sum += A.get(i).replaceAll("-","").length();
		}
		return 1 - Distance.multiDistance(A, B)/ (double)( (A.size()-1) *sum );
	}

	public static double sequenceSimilarityScore(String []  A){
		ArrayList<String> temp = new ArrayList<String>(Arrays.asList(A)); 
		return sequenceSimilarityScore(temp);
	}

	public static double sequenceSimilarityScore(List<String>  A){
		double score = 0;
		int count = 0;
		int tempScore;
		int length;
		for(int i =0; i<A.size()-1; ++i){
			for(int j =i+1; j<A.size(); ++j){
				String seq1 = A.get(i);
				String seq2 = A.get(j);
				tempScore = 0;
				length = 0;
				count++;
				for(int k = 0; k<seq1.length(); ++k){
					if (seq1.charAt(k) == '.' || seq2.charAt(k) == '.'){
						continue;
					}
					else if(seq1.charAt(k) == seq2.charAt(k)){
						tempScore++;
						length++;
					}else{
						length++;
					}
				}
				score += tempScore / (double)length;
			}
		}
		return score /count;
	}



}




