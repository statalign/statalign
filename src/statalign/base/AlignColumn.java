package statalign.base;

/**
 * 
 * An AlignColumn is an object that represents an observed or a hidden character
 * and stores a double array storing the Felsenstein likelihoods.
 * 
 * AlignColumns are double chained to form an observed or hidden sequence.
 * There are also links to their parent. If the character that AlignColumn represents 
 * is aligned to an ancestral characters, then the alignColumn is not <tt>orphan</tt>
 * otherwise it is, and then its  <tt>parent</tt> represents the first ancestral character
 * to the right in the alignment column.
 * 
 * @author miklos, novak
 *
 */
public class AlignColumn {

	/**
	 * The owner tells the tree vertex where this AlignColumn belongs to. 	
	 */
	public Vertex owner;

	/**
	 * This reference points to the parent of this AlignColumn
	 */
	public AlignColumn parent;
	/**
	 * <p>If true, then this AlignColumn is inserted, and its <tt>parent</tt> is the
	 * first ancestral character to the right in the alignment
	 * 
	 * <p>If it is false, then this AlignColumn is aligned to its parent 
	 */
	public boolean orphan;					/* if false, column has a true parent
														 if true, parent shows parent vertex's
														 next column in alignment */
    boolean selected;     	  /* if true, it is selected as a homologous column that is not changed
    								during changing topology
     						*/

    boolean emptyWindow;   /* if true, then it is the end of an empty window */

    /**
     * This reference points to the previous AlignColumn in the alignment
     */
	public AlignColumn prev;				// previous column in alignment
	/**
	 * This reference points to the next AlignColumn in the alignment
	 */
	public AlignColumn next;				// next column in alignment

	/**
	 * This reference points to the left child of this AlignColumn
	 */
	public AlignColumn left;				// left child in alignment
	/**
	 * This reference points to the right child of this Aligncolumn
	 */
	public AlignColumn right;			// right child in alignment

	/**
	 * This double array contains the Felsenstein Likelihoods
	 */
	public double seq[];						// Felsenstein likelihoods of the column
	/**
	 * This double array contains the upper Likelihoods
	 */
	public double upp[];						// upper likelihoods of the column

	/**
	 * It constructs a new AlignColumn. Sets only the owner, other fields are filled in outside of the
	 * constructor
	 * @param owner The vertex where this AlignColumn belongs to.
	 */
	public AlignColumn(Vertex owner){
		this.owner = owner;
		parent = null;
		orphan = true;
	}
	
	public AlignColumn clone() {
		AlignColumn c = new AlignColumn(owner);
		c.parent = parent;
		c.orphan = orphan;
		c.left = left;
		c.right = right;
		c.selected = selected;
		c.emptyWindow = emptyWindow;
		if (seq==null) c.seq = null;
		else 		   c.seq = seq.clone();
		return c;
	}
	/**
	 * It creates a new AlignColumn, chains it to the next column (namely, it is used to
	 * generate a new ancestral sequence built in a traceback phase of a dynamic programming).
	 * the owner is set to the same vertex as the owner of the next AlignColumn.
	 * 
	 * @param next The next AlignColumn in the sequence that is being generated
	 * @param newSeq If true, then the constructor allocates memory for
	 * the seq[] array that contains the Felsenstein likelihoods. If false, then 
	 * the seq[] array of other AlignColumn is associated to this Aligncolumn and reused.
	 * For speeding-up purposes.		
	 */
	public AlignColumn(AlignColumn next, boolean newSeq) {
		owner = next.owner;
		parent = null;
		orphan = true;
		this.next = next;
		next.prev = this;
		if(newSeq) {
			seq = new double[owner.owner.substitutionModel.e.length];
			upp = new double[owner.owner.substitutionModel.e.length];
		}
	}	
	public AlignColumn(AlignColumn next) {
		owner = next.owner;
		parent = null;
		left = null;
		right = null;
		orphan = true;
		this.next = next;
		prev = next.prev;
		if (prev!=null) prev.next = this;
		next.prev = this;
		seq = new double[owner.owner.substitutionModel.e.length];
		upp = new double[owner.owner.substitutionModel.e.length];		
		owner.length++;
	}

	/** 
	 *  Remove all non-parent references to a column, <code>p</code>.
	 *  Any columns that called this column 'parent' will have a null
	 *  parent after this function is called, and must be adopted
	 *  (unless they belong to the root Vertex).
	 *  @param p The column to be deleted.
	 */
	static void delete(AlignColumn p) {
		
//		if (Utils.DEBUG) {
//			if (p.left!=null && p.right!=null) {
//				throw new IllegalArgumentException("A column cannot be deleted if it has two descendants.");
//			}
//		}
		
		if (Utils.DEBUG) System.out.println("Deleting column "+p.toString()+" at Vertex "+p.owner.index);
		if (!p.orphan) {
			if (p==p.parent.left)  p.parent.left = null;
			else  				   p.parent.right = null;			
		}

		// Now skip column p in the parent sequence
		if (p.prev != null) p.prev.next = p.next;
		if (p.next != null) p.next.prev = p.prev;     
		
		if (p.left != null)  { p.left.orphan = true; p.left.parent = null; }
		if (p.right != null) { p.right.orphan = true; p.right.parent = null; }
		// NB the reason for making the new orphan's parent field null
		// is to allow the memory associated with p to be freed
		// as soon as possible, rather than waiting for p.left and p.right
		// to be adopted by another parent.
		
		p.owner.length--;
	}
	
	/**
	 * Updates the parent/child pointers with a specified parent.
	 * @param p The parent column. 
	 * @param isLeft Whether <code>this</code> is to be the left descendant of <code>p</code>.
	 */
	void updateParent(AlignColumn p, boolean isLeft) {		
		
		// Remove this node from the list of children of its old parent
		if (p != parent) {
			if (parent!=null) {
				if (this==parent.left) parent.left = null;
				else if (this==parent.right) parent.right = null;
				else ; // the original parent had already taken a new child
			}
			parent = p;
		}
		
		//if (orphan) parent = p.next; 
		//else {		
		if (!orphan) {			
			if (isLeft) p.left = this;
			else 	  p.right = this;
		}
	}
	
	void setWinFirst(AlignColumn prev) {
		owner.winFirst = this;
		this.prev = prev;
		if(prev != null)
			prev.next = this;
		else
			owner.first = this;
	}

	void saveLeft(AlignColumn copy) {
		if(copy.left != null) {
			left = copy.left;
		}
		selected = copy.selected;
		emptyWindow = copy.emptyWindow;
	}

	void saveRight(AlignColumn copy) {
		if(copy.right != null) {
			right = copy.right;
		}
		selected = copy.selected;
		emptyWindow = copy.emptyWindow;
	}

	void saveBoth(AlignColumn copy) {
		saveLeft(copy);
		saveRight(copy);
		seq = copy.seq;
		upp = copy.upp;
	}

	void saveParent(AlignColumn copy) {
		if(copy.parent != null) {
			parent = copy.parent;
			orphan = copy.orphan;
			if(!orphan) {
				if(parent.left == copy)
					parent.left = this;
				else
					parent.right = this;
			}
		}
	}

	AlignColumn brother() {
		if(parent != null && !orphan) {
			return parent.left == this ? parent.right : parent.left;
		}
		return null;
	}

	char mostLikely(){
	//	double max = 0.0;
	//	char character = '*';
	//	for(int i = 0; i < seq.length; i++){
	//		if(seq[i] > max){
	//			character =  owner.owner.substitutionModel.alphabet[i];
	//			max = seq[i];
				//			s[1] += owner.substitutionModel.alphabet[i];
				//found = true;
	//		}
	//	}

	//	return character;
		return owner.owner.substitutionModel.mostLikely(seq);
	}
/**
    This function tells if a character is homolgous in all selected vertices
 */
 boolean isHomologous(){
	boolean x = true;
	if(owner.left == null){
	    return true;
	}
	if(owner.left.selected && left != null){
	    x = x && left.isHomologous();
	}
	if(owner.left.selected && left == null){
	    return false;
	}
	if(owner.right.selected && right != null){
	    x = x && right.isHomologous();
	}
	if(owner.right.selected && right == null){
	    return false;
	}
	
	return x;
 }
 
 /**
    This function tells if it is an empty window in all selected vertices
    Works only for selected AlignColumns, since we dom't check if left and right exist
  */
  boolean isEmptyWindow(){ 
	boolean x = (prev == null || prev.selected);
	if(owner.left != null && owner.left.selected){
	    x = x && left.isEmptyWindow();
	}
	if(owner.right != null && owner.right.selected){
	    x = x && right.isEmptyWindow();
	}
	return x;
 }

 void selectDown(){
	selected = true;
	if(owner.left != null && owner.left.selected){
	    left.selectDown();
	}
	if(owner.right != null && owner.right.selected){
	    right.selectDown();
	}
 }

 void setEmptyWindowDown(boolean value){
	emptyWindow = value;
	if(owner.left != null && owner.left.selected){
	    left.setEmptyWindowDown(value);
	}
	if(owner.right != null && owner.right.selected){
	    right.setEmptyWindowDown(value);
	}

 }

 /**
    This function finds the beginning of the window in fast topology changing
  */
  AlignColumn findWindowStart(){
	AlignColumn a = this;
	int length = 0;
	//	System.out.println("a: "+a+" selected: "+a.prev.selected);
	while(a.prev != null && !(a.prev.selected && a.prev.emptyWindow)){
	    a = a.prev;
	    length ++;
	}
	owner.winLength = length;
	return a;
 }

}
