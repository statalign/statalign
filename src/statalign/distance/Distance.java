package statalign.distance;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;



public class Distance {
	
	

	public static void main(String[] args) {

		ArrayList<String> A1 = new ArrayList<String>();
		ArrayList<String> A2 = new ArrayList<String>();
		ArrayList<String> amapseq = new ArrayList<String>();
		

		A1.add("GGGUGCUUGAAGCUGUCUGCUUUAAGUGCUUGCA----UCAGGCUGAGAGUAGGCAGAGAAAAGCCCCGUA------------------UCA-----A----------------UGUUAAUCAAUACGAGGC-CCUCUGUAAUG-CACGACAACAUUACGGU-AGCCUUUUACC-CGCCGAAA-GGCAA------GGAGGCUGAAGAUG");
		A1.add("GGGUGCUUGAGACUGUUUGUCUCAGG------UAUUUA----CCAAAAGGCAGACAGAGAAAAGCCCCACC------------------UGACUAUA------------------AAUCAAAAGUGCAUUGC-ACCCAUUAUUG-AUCU-CUUCAAUAACGA-AGCUA-UCCCC-UACAGUAU-UUCA------GAACGUCCAACCAUG");
		A1.add("GGGUGCUUGAGGCUGUCUGCCUCGGG------CAUGCC---ACUGUAAGGCAGACAGAGAAAAGCCCCAGUUAACAUUACGCGUCCUGC--------AAGA----------CGCUUAACAUUAAUCUGAGGC-CCAAUCUAUGU-CUCA-CAAAUGUAGGUU-AGCCUCUUACG-UGCCGAAA-GGCAAGGAGAAGCAGGCUAUG-AAG");
		A1.add("GGGUGCUUGAGACUGUUUGUCUCAGG------UAUUCA----CCGAAAGGCAGACAGAGAAAAGCCCCACC------------------UGACU---------------------AUAAAUCAAAGUGAGGCUA--CCCU-AUGCCUGAACACCAUAAGG-UUAGCCUCUUACUCGUUGGAAAUCAACAC----AGGGGGCUGGGAAUG");
		A1.add("GGGUGCUUGAGGCUGUCUGCCUCGGG------CAUGCC---ACCGUAAGGCAGACAGAGAAAAGCCCCAGU------------------UAACAUUACGCGUCCUGCAAGACGCCUAACAUUAAUCUGAGGC-CAAUUU-CAUG-CUAGACA-CAUGUAGGUUAGCCUCUUACG-CGCCGAAA-GGCAAG----GAGAAGCAGCU-AUG");

		A2.add("GGGUGCUUGAGGCUGUCUGCCUCGGGCA--UGCCAC---UGUAAGGCAGACAGAGAAAAGCCCCAGUU-AACAUUACGCGUCCUGCAAGACGCUUAACAUUAA-UCUGAGGCCCAAUCU-AUGUCUCA-CAAA---UGU---AGGUUAGCCUCUUACGUGCCGAAAGGCAA----GGAGAAGCAGGCU-AUGAAG");
		A2.add("GGGUGCUUGAAGCUGUCUGCUUUAAGUGCUUGCAUCAGGCUGAGAGUAGGCAGAGAAAAGCCCCGUAUCAAUGUUA----------------------AUCAA-UACGAGGCCCUCUGUAAUGCACGA-CAAC---AUU---ACGGUAGCCUUUUACCCGCCGAAAGGCAA----GGAG------GCUGAAGAUG");
		A2.add("GGGUGCUUGAGACUGUUUGUCUCAGGUAUU---UAC---CAAAAGGCAGACAGAGAAAAGCCCCACCUGACUAUAA----------------------AUCAAAAGUGCAUUGCAC-CC-AUUAUUGAUCUCUUCAAUAACGAAGCUAUCCCCU-----ACAGUAUUUCAG----AACG------UCCAACCAUG");
		A2.add("GGGUGCUUGAGGCUGUCUGCCUCGGGCA--UGCCAC---CGUAAGGCAGACAGAGAAAAGCCCCAGUU-AACAUUACGCGUCCUGCAAGACGCCUAACAUUAA-UCUGAGGCCAAUUUC-AUGCUAGA-CACA---UGU---AGGUUAGCCUCUUACGCGCCGAAAGGCAA----GGAGAAGCA-GCU-AUG---");
		A2.add("GGGUGCUUGAGACUGUUUGUCUCAGGUAUU---CAC---CGAAAGGCAGACAGAGAAAAGCCCCACCUGACUAUAA----------------------AUCAA-AGUGAGGCU-ACCCU-AUGCCUGAACACC---AUA---AGGUUAGCCUCUUACUCGUUGGAAAUCAACACAGGGG------GCUGGGAAUG");
		
		
		
		amapseq.add("GGGUGCUUGAAGCUGUCUGCUUUAAGUGCUUGCAUCAGGCUGAGAGUA----------------------------------------------------------------------------------------------------------------------GGCAGAGAAAAGCCCCGUAUCAAUGUUAAUCAAUACGAGGCCCUCUGUAAUGCACGACAACAUUACGGUAGCCUUUUACCCGCCGAAAGGCAAGGAGGCUGAAGAUG------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
		amapseq.add("GGGUGCUUGAGAC-----------------------------------UGUUUGUCUCAGGUAUUUACCAAAAGGCA-----------------------------------------------------------------------------------------GACAGAGAAAAGCCCC-------------------------------------------------------------------------------------------ACCUGACUAUAAA-------------------------------------------------------------------------------------------------------------UCAAAAGUGCAUUGCACCCAUUAUUGAUCUCUUCAAUAACGAAGCUAUCCCCUACAGUAUUUCAGAACGUCCAACCAUG---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
		amapseq.add("GGGUGCUUGAGGC---------------------------------------------------------------------------------------------------------------------------UGUCUGCCUCGGGCAUGCCACUGUAAGGCAGACAGAGAAAAGCCCC--------------------------------------------------------------------------------------------------------AGUUAACAUUACGCGUCCUGCAAGACGCUUAACAUUAAUCUGAGGCCCAAUCUAUGUCUCACAAAUGUAGGUUAGCCUCUUACGUGCCGAAAGGCAAGGAGAAGCAGGC-------------------------------------------------------------------------------------------------------------------------------------UAUGAAG--------------------------------------------------------------------------------------------------------------------------------------------------------");
		amapseq.add("GGGUGCUUGAGAC----------------------------------------------------------------UGUUUGUCUCAGGUAUUCACCGAAAGGCA------------------------------------------------------------GACAGAGAAAAGCCCC----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ACCUGACUAUAAAUCAAAGUGAGGCUACCCUAUGCCUGAACACCAUAAGGUUAG-------CCUCUUACUCGUUGGAAAUCAACACAGGGGGCUGGGAAUG----------------------------------------------------------------------------------------------------------------");
		amapseq.add("GGGUGCUUGAGGC---------------------------------------------------------------------------------------------UGUCUGCCUCGGGCAUGCCACCGUAAGGCA------------------------------GACAGAGAAAAGCCCC---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AGUUAACAUUACGCGUCCUGCAAGACGCCUAACAUUAAUCUGAGGCCAAUUUCAUGCUAGACACAUGUAGGUUAGCCUCUUACGCGCCGAAAGGCAAGGAGAAGCAGCUAUG");
		
		
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

		//A1[0] = "ATGCTAC";
		//A1[1] = "A--CTAC";
		//A1[2] = "A-----C";

		//A2[0] = "ATG-C-TAC";
		//A2[1] = "A---CTA-C";
		//A2[2] = "A------C-";

		System.out.println("The distance between the reference alignment and the MPD from Test1: "+Distance.multiDistance(A1,A2));
		//System.out.println(Distance.distance(A1.toArray(),A2.toArray()));
		System.out.println("AMA between those alignments: "+Distance.AMA(A1,A2));
		System.out.println("-------------------------------------------");
		System.out.println("The distance between the reference alignment and the AMAP from Test1: "+Distance.multiDistance(A1,amapseq));
		//System.out.println(Distance.distance(A1.toArray(),A2.toArray()));
		System.out.println("AMA between those alignments: "+Distance.AMA(A1,amapseq));
		System.out.println("-------------------------------------------");
		System.out.println("The distance between the AMAP and the MPD from Test1: "+Distance.multiDistance(amapseq,A2));
		//System.out.println(Distance.distance(A1.toArray(),A2.toArray()));
		System.out.println("AMA between those alignments: "+Distance.AMA(amapseq,A2));
		System.out.println("-------------------------------------------");
		//System.out.println("Distance between the reference alignment and the MPD from Test17: "+Distance.multiDistance(B1,B2));
		//System.out.println(Distance.distance(A1.toArray(),A2.toArray()));
		//System.out.println("AMA between those alignments: "+Distance.AMA(B1,B2));






	}
	
	/**
	 * Given a set of sequences and a similarity score, returns the corresponding distance.
	 * @param sequences
	 * @param sim
	 * @return
	 */
	public static double similarityToDifference(ArrayList<String> sequences, double sim)
	{
		int lengthSum = 0;
		for(int i = 0 ; i < sequences.size() ; i++)
		{
			lengthSum += sequences.get(i).replaceAll("-", "").length();
		}
		int ksub1 = sequences.size() - 1;
		
		return (1-sim)*ksub1*lengthSum;
	}

	public static int multiDistance(ArrayList<String>  A, ArrayList<String>  B){
		sortSeq(A,B);
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

	public static double AMA(ArrayList<String>  A, ArrayList<String>  B){
		int sum = 0;
		for(int i = 0; i<A.size(); ++i){
			sum += A.get(i).replaceAll("-","").length();
		}
		return 1 - Distance.multiDistance(A, B)/ (double)( (A.size()-1) *sum );
	}

	
}




