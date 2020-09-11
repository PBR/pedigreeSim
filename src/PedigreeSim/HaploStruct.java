/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Locale;
import java.util.TreeSet;

/**
 * HaploStruct describes how a chromosome is composed of segments of the
 * original founder alleles. This information is stored in two ArrayLists:
 * - founder is a sequence of integers, each identifying a founder
 *   allele; for each segment there is one integer element;
 * - recombpos lists the positions of the recombinations as a sequence
 *   of doubles. For each segment (each element of founder) there is also
 *   a position in recombPos: this is the position where the segment STARTS.
 *   The first recombPos entry is the chromosome's headPos. The end of the
 *   final segment is not listed in recombPos, but is equal to the
 *   chromosome's tailPos.
 * A HaploStruct without recombinations has only one segment, i.e. one
 * entry in the founder and recombPos ArrayLists.
 * @author Roeland Voorrips
 */
public class HaploStruct implements Cloneable {
    private ArrayList<Integer> founder;
        //founder allele of each segment
    private ArrayList<Double> recombPos;
        //recombination positions, start position of each segment, first is chrom.headPos
    private Chromosome chrom;

    //protected NumberFormat = null;

    public HaploStruct(Chromosome chromosome, int startFounderAllele) 
            throws Exception {
        if (chromosome==null) {
            throw new Exception("Chromosome==null in HaploStruct constructor");
        }
        chrom = chromosome;
        recombPos = new ArrayList<Double>();
        founder = new ArrayList<Integer>();
        recombPos.add(chrom.getHeadPos());
        founder.add(startFounderAllele);
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        try {
            HaploStruct copy = (HaploStruct) super.clone();
            copy.chrom = chrom; //is never modified, no cloning needed
            //make a deep copy of the arraylists:
            copy.recombPos = new ArrayList<Double>();
            copy.founder = new ArrayList<Integer>();
            for (int i=0; i<founder.size(); i++) {
                //conversion to doublevalue and intvalue to force a new
                //Double or Integer object
                copy.recombPos.add(recombPos.get(i));
                copy.founder.add(founder.get(i));
            }
            return copy;
        }
        catch( CloneNotSupportedException e ) {
            return null;
        }
    } //clone


    public ArrayList<Integer> getFounder() {
        return founder;
    }

    public ArrayList<Double> getRecombPos() {
        return recombPos;
    }

    public Chromosome getChrom() {
        return chrom;
    }

    /**
     * toString: writes a Haplostruct as an sequence of founder alleles
     * (integers) alternating with recombination points (doubles)
     * separated by tabs
     * @return 
     */
    @Override
    public String toString() {
        DecimalFormat fix = new DecimalFormat("#0.0000",new DecimalFormatSymbols(Locale.US));
        String s = founder.get(0).toString();
        for (int i=1; i<founder.size(); i++) {
            s += "\t" + fix.format(recombPos.get(i)) +
                 "\t" + founder.get(i);
        }
        return s;
    }

    /**
     * findSegment :
     * finds the last segment of which the start position is <= the argument
     * @param position
     * @return
     */
    public int findSegment(double position) {
        int seg = recombPos.size()-1;
        while (seg>0 && recombPos.get(seg)>position)
            seg--;
        return seg;
    }

    /**
     * addSegment adds a new segment at the end
     * @param position (in Morgan) where the new segment starts
     * @param founderAllele is the new segment
     */
    public void addSegment(double position, int founderAllele) {
        recombPos.add(position);
        founder.add(founderAllele);
    }

    /**
     * insertSegment inserts a new segment starting at position
     * @param position (in Morgan) where the new segment starts
     * @param founderAllele is the new segment
     */
    public void insertSegment(double position, int founderAllele) {
        int seg = findSegment(position);
        recombPos.add(seg+1, position);
        founder.add(seg+1, founderAllele);
    }

    public int getFounderAt(double position) throws Exception {
        if (position<chrom.getHeadPos() ||
            position>chrom.getTailPos()) {
            throw new Exception("getFounderAt: invalid position");
        } else {
            return founder.get(findSegment(position));
        }
    }

    public int getFounderAtCentromere() {
        int f=0;
        try {
            f = getFounderAt(chrom.getCentromerePos());
        } catch (Exception ex) {
            System.out.println("getFounderAtCentromere: "+ex.getMessage());
        }
        return f;
    }

    public int getFounderAtHead() {
        int f=0;
        try {
            f = getFounderAt(chrom.getHeadPos());
        } catch (Exception ex) {
            System.out.println("getFounderAtHead: "+ex.getMessage());
        }
        return f;
    }

    public int getFounderAtTail() {
        int f=0;
        try {
            f = getFounderAt(chrom.getTailPos());
        } catch (Exception ex) {
            System.out.println("getFounderAtTail: "+ex.getMessage());
        }
        return f;
    }

    public int segmentCount() {
        return founder.size();
    }

    /**
     * internal function used by recombine()
     * @param hs
     * @param recombPos
     * @return
     * @throws Exception 
     */
    private boolean recombinationInitialization(HaploStruct hs, double recombPos)
            throws Exception {
        if (hs==null || !this.getChrom().equals(hs.chrom)) {
            throw new Exception("recombine: parameter hs invalid");
        }
        if (recombPos<chrom.getHeadPos() || recombPos>chrom.getTailPos()) {
            throw new Exception("recombine: parameter recombPos invalid");
        }
        return recombPos>chrom.getHeadPos() && recombPos<chrom.getTailPos();
            //recombination exactly at head or tail will be ignored)
    }

    /**
     * recombine_keepHead produces the result of a recombination at position
     * recombPos between HaploStructs this and hs.
     * The recombined Haplostructs are reassigned to this and hs such that
     * the heads of the chromosomes stay where they were.
     * @param hs must belong to the same chromosome as this
     * @param recombPos a recombination position outside the chromosome head
     * and tail throws an Exception; recombinations exactly at the head and tail
     * are ignored.
     * @throws Exception
     */
    public void recombine_keepHead(HaploStruct hs, double recombPos) throws Exception {
        if (recombinationInitialization(hs, recombPos)) {
            //rombination exactly at head or tail is ignored)

            //find the segments of each where the chiasma is located:
            int thisseg = this.findSegment(recombPos);
            int hsseg = hs.findSegment(recombPos);
            //copy all segments starting from thisseg from this to tmp,
            //and remove all starting at thisseg+1 from this:
            ArrayList<Double> tmprecpos = new ArrayList<Double>();
            ArrayList<Integer> tmpfounder = new ArrayList<Integer>();
            tmprecpos.add(recombPos);
            tmpfounder.add(this.founder.get(thisseg));
            while (this.founder.size()>thisseg+1) {
                tmprecpos.add(this.recombPos.get(thisseg+1));
                this.recombPos.remove(thisseg+1);
                tmpfounder.add(this.founder.get(thisseg+1));
                this.founder.remove(thisseg+1);
            }
            //copy all segments starting from hsseg from hs to this,
            //and remove all starting at hsseg+1 from hs:
            this.recombPos.add(recombPos);
            this.founder.add(hs.founder.get(hsseg));
            while (hs.founder.size()>hsseg+1) {
                this.recombPos.add(hs.recombPos.get(hsseg+1));
                hs.recombPos.remove(hsseg+1);
                this.founder.add(hs.founder.get(hsseg+1));
                hs.founder.remove(hsseg+1);
            }
            //copy all segments from tmp to hs:
            for (int s=0; s<tmpfounder.size(); s++) {
                hs.recombPos.add(tmprecpos.get(s));
                hs.founder.add(tmpfounder.get(s));
            }
        } //recombPos between (not at) chrom head and tail
    } //recombine_keepHead

    /**
     * recombine_keepTail produces the result of a recombination at position
     * recombPos between HaploStructs this and hs.
     * The recombined Haplostructs are reassigned to this and hs such that
     * the ends of the chromosomes stay where they were.
     * @param hs must belong to the same chromosome as this
     * @param recombPos a recombination position outside the chromosome head
     * and tail throws an Exception; recombinations exactly at the head and tail
     * are ignored.
     * @throws Exception
     */
    public void recombine_keepTail(HaploStruct hs, double recombPos) throws Exception {
        if (recombinationInitialization(hs, recombPos)) {
            recombine_keepHead(hs, recombPos);
            /* now the chromosome ends are switched and must be switched back:
             */
            ArrayList<Double> tmprecpos; // = new ArrayList<Double>();
            ArrayList<Integer> tmpfounder; // = new ArrayList<Integer>();
            tmprecpos = this.recombPos;
            tmpfounder = this.founder;
            this.recombPos = hs.recombPos;
            this.founder = hs.founder;
            hs.recombPos = tmprecpos;
            hs.founder = tmpfounder;
        }
    } //recombine_keepTail

    /**
     * recombine produces the result of a recombination at position recombPos
     * between HaploStructs this and hs.
     * The recombined Haplostructs are reassigned to this and hs such that
     * the centromeres stay where they were.
     * @param hs must belong to the same chromosome as this
     * @param recombPos a recombination position outside the chromosome head
     * and tail throws an Exception; recombinations exactly at the head and tail
     * are ignored.
     * @throws Exception
     */
    public void recombine(HaploStruct hs, double recombPos) throws Exception {
        recombine_keepHead(hs, recombPos);
        /* if the recombination took place left of the centromere, now the
         * centromeres are switched and must be switched back:
         */
        if (recombinationInitialization(hs, recombPos) &&
                recombPos<chrom.getCentromerePos()) {
            ArrayList<Double> tmprecpos; // = new ArrayList<Double>();
            ArrayList<Integer> tmpfounder; // = new ArrayList<Integer>();
            tmprecpos = this.recombPos;
            tmpfounder = this.founder;
            this.recombPos = hs.recombPos;
            this.founder = hs.founder;
            hs.recombPos = tmprecpos;
            hs.founder = tmpfounder;
        }
    } //recombine

    /**
     * getHaplotypeString :
     * lists the actual alleles at all loci along this HaploStruct
     * @return 
     */
    public String getHaplotypeString() {
        ArrayList<Locus> locus = chrom.getLocus();
        if (locus==null || locus.isEmpty()) {
            return "";
        } else {
            String s="";
            for (Locus loc: locus) {
                if (!s.isEmpty()) s += "\t";
                int f = founder.get(findSegment(loc.getPosition()));
                s += loc.getAlleleName(f);
            }
            return s;
        }
    } //getHaplotypeString

    /**
     *
     * @param position
     * @return true if founderallele at position is different from that
     * at head of chromosome
     * @throws Exception
     */
    public boolean isRecombinantHead(double position) throws Exception {
        return getFounderAt(position) != getFounderAt(chrom.getHeadPos());
    }

    /**
     *
     * @param position
     * @return true if founderallele at position is different from that
     * at tail of chromosome
     * @throws Exception
     */
    public boolean isRecombinantTail(double position) throws Exception {
        return getFounderAt(position) != getFounderAt(chrom.getTailPos());
    }

    public TreeSet<Integer> founderAlleles() {
        TreeSet<Integer> alleles = new TreeSet<Integer>();
        for (Integer f : founder) {
            alleles.add(f);
        }
        return alleles;
    } //founderAlleles

    /**
     * @return the number of different founder alleles represented 
     * in the HaploStruct
     */
    public int founderCount() {

        int[] present = new int[chrom.getPopdata().founderAlleleCount];
        for (int i=0; i<chrom.getPopdata().founderAlleleCount; i++) {
            present[i] = 0;
        }
        for (Integer founder1 : founder) {
            present[founder1] = 1;
        }
        int sum = 0;
        for (int i=0; i<chrom.getPopdata().founderAlleleCount; i++) {
            sum += present[i];
        }
        return sum;
    } //founderCount

}
