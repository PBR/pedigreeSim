/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;
import java.util.TreeSet;

/**
 * Genotype is the base class for Individuals (which can have a name and parents) 
 * and for Gametes (which can do fertilization)
 * @author Roeland Voorrips
 */
public class Genotype {

     /**
     * Class ChromConfig describes for one chromosome in one meiosis
     * - quad: how many Quadrivalents are formed
     * - chromseq: the order of the homologs:
     *     from the start, quad groups of 4 homologs form the quadrivalents;
     *     with cross-type quadrivalents the 1st and 4th, and the 2nd and 3rd
     *     are paired at the head/north ends while the 1st and 2nd, and
     *     the 3rd and 4th are paired at the tail/south ends.
     *     After that the (ploidy-4*quad)/2 pairs form the bivalents.
     * Used by the classes Individual and Gamete:
     *   Individual: to get the meiotic configurations in doMeiosis,
     *               to transmit this to the Gametes formed,
     *               and to store the configurations of the two meioses
     *               from which it originated
     *   Gamete: to store the configurations of the meiosis
     *               from which it originated, and
     *               to transmit that to the new Individual in fertilisation
     */
    static class ChromConfig {
        int quad; //the number of quadrivalents
        int[] chromseq;

        protected ChromConfig(int ploidy) {
            chromseq = new int[ploidy];
        }

        protected ChromConfig(ChromConfig chrconf) {
            quad = chrconf.quad;
            chromseq = new int[chrconf.chromseq.length];
            System.arraycopy(chrconf.chromseq, 0, chromseq, 0, chromseq.length);
        }
    }

    protected PopulationData popdata;
    /**
     * haplostruct has two indices:
     * first index 0..popdata.chromosomeCount-1,
     * second index 0..ploidy-1
     */
    protected HaploStruct[][] haplostruct;

    //protected fields for faster access, must be set by descendant classes:
    protected Tools tools;
    protected Random rand;

    /**
     * NOTE: if this constructor is used, haplostruct must be set later !!!
     * @param popdata
     * @throws Exception
     */
    protected Genotype(PopulationData popdata) throws Exception {
        if (popdata==null) {
            throw new Exception("Genotype constructor: popdata==null");
        }
        this.popdata = popdata;
        if (popdata.tools==null) {
            throw new Exception("Genotype constructor: popdata.tools==null");
        }
        this.tools = popdata.tools;
        if (popdata.tools.rand==null) {
            throw new Exception("Genotype constructor: popdata.tools.rand==null");
        }
        this.rand = popdata.tools.rand;
    }

    /** This constructor checks if haplostruct has the correct dimensions
     * and if so, assigns haplostruct to this.haplostruct; 
     * else throws an Exception
     * @param ploidy The expected ploidy of this Genotype (different for Individual
     * and for Gamete)
     * @param haplostruct
     * @param popdata
     * @throws Exception
     */
    protected Genotype(int ploidy, HaploStruct[][] haplostruct,
            PopulationData popdata)
            throws Exception {
        this(popdata);
        if (!checkHaploStruct(ploidy,haplostruct)) {
            throw new Exception("haplostruct invalid in Genotype constructor");
        }
        this.haplostruct = haplostruct;
    }

    /**
     * For use by constructors of descendant classes
     *
     * @param ploidy the expected ploidy of the Genotype (=popdata.ploidy
     *        for an Individual, half of that for a Gamete
     * @param haplostruct first index 0..popdata.chromosomeCount-1,
     *        second index 0..ploidy-1
     * @return
     */
    protected boolean checkHaploStruct(int ploidy, HaploStruct[][] haplostruct) {
        boolean ok = (haplostruct != null &&
                      haplostruct.length==popdata.chromCount());
        int c=0;
        while (ok && c<popdata.chromCount()) {
            ok = haplostruct[c]!=null && haplostruct[c].length==ploidy;
            c++;
        }
        return ok;
    }

    public HaploStruct[][] getAllHaploStruct() {
        return haplostruct;
    }

    public HaploStruct getHaploStruct(int chromnr, int homolog) {
        return haplostruct[chromnr][homolog];
    }
    
    public void clearAllHaplostruct() {
        haplostruct = null;
    }

    public PopulationData getPopdata() {
        return popdata;
    }

    /**
     * getLocusAllele returns a String array with the getPloidy() names
     * of the alleles in this Genotype at the specified chromosome and locus
     * @param c the chromosome number
     * @param loc
     * @return
     * @throws Exception 
     */
    public String[] getLocusAllele(int c, int loc) throws Exception {
        String[] alleles = new String[getPloidy()];
        for (int i=0; i<alleles.length; i++) {
            Locus locus =  popdata.getChrom(c).getLocus().get(loc);
            int founder = haplostruct[c][i].getFounderAt(locus.getPosition());
            alleles[i] = locus.getAlleleNames()[founder];
        }
        return alleles;
    }

    /**
     * getLocusAllele returns a String array with the getPloidy() names
     * of the alleles in this Genotype at the specified chromosome and locus
     * @param c the chromosome number
     * @param pos position
     * @return
     * @throws Exception 
     */
    public int[] getFounderAlleles(int c, double pos) throws Exception {
        int[] founderAlleles = new int[getPloidy()];
        for (int i=0; i<founderAlleles.length; i++) {
            founderAlleles[i] = haplostruct[c][i].getFounderAt(pos);
        }
        return founderAlleles;
    }

    /**
     * @param c the chromosome number
     * @param founder if true, homozygosity for founder alleles will be
     * checked, else homozygosity for actual (observed) alleles.
     * @return a boolean[] with as  index the locus number on the chromosome.
     * The values are true if homozygous (same locus allele on all homologous
     * chromosomes) and false if heterozygous.
     */
    public boolean[] homozygous(int c, boolean founder) {
        boolean[] result;
        Chromosome chrom = popdata.getChrom(c);
        if (chrom==null || chrom.getLocus()==null ||
                chrom.getLocus().isEmpty()) {
            result = null;
        }
        else {
            result = new boolean[chrom.getLocus().size()];
            int genplo = getPloidy(); //ploidy of this Genotype
            Locus locus;
            double locPos;
            int firstFounder;
            String first = null;
            for (int loc=0; loc<result.length; loc++) {
                try {
                    locus = chrom.getLocus().get(loc);
                    locPos = locus.getPosition();
                    firstFounder = haplostruct[c][0].getFounderAt(locPos);
                    if (!founder) first = locus.getAlleleName(firstFounder);
                    result[loc] = true;
                    int h = 1;
                    while (result[loc] && h<genplo) {
                        if (founder) {
                            result[loc] = firstFounder == 
                                haplostruct[c][h].getFounderAt(locPos);
                        } else {
                            result[loc] = first.equals(locus.getAlleleName(
                                haplostruct[c][h].getFounderAt(locPos)));
                        }    
                        h++;
                    }
                }  catch (Exception ex) {/* cannot occur */ }
            }
        }
        return result;
    } //homozygous

    /**
     * @param c the chromosome number
     * @return an int[][] with as first index the locus number on the chromosome
     * and second index the founder alleles in the population
     * The values are the dosages of each founder allele at each locus.
     */
    public int[][] getDosages(int c) {
        int[][] result;
        Chromosome chrom = popdata.getChrom(c);
        if (chrom==null || chrom.getLocus()==null ||
                chrom.getLocus().isEmpty() ||
                popdata.founderAlleleCount==0) {
            result = null;
        }
        else {
            result = new int[chrom.getLocus().size()][popdata.founderAlleleCount];
            //initialized containing zeroes
            int genplo = getPloidy(); //ploidy of this Genotype
            Locus locus;
            double locPos;
            for (int loc=0; loc<chrom.getLocus().size(); loc++) {
                try {
                    locus = chrom.getLocus().get(loc);
                    locPos = locus.getPosition();
                    for (int h=0; h<genplo; h++) {
                        int f = haplostruct[c][h].getFounderAt(locPos);
                        result[loc][f] += 1;
                    }
                }  catch (Exception ex) {/* cannot occur */ }
            }
        }
        return result;
    } //getDosage

    public int getPloidy() {
        if (haplostruct==null || haplostruct.length==0 ||
            haplostruct[0]==null) {
            return 0;
        }
        else {
            return haplostruct[0].length;
        }
    }

    public TreeSet<Integer> founderAlleles() {
        TreeSet<Integer> alleles = new TreeSet<Integer>();
        for (int c =0; c<popdata.chromCount(); c++) {
            for (int h=0; h<getPloidy(); h++) {
                alleles.addAll(haplostruct[c][h].founderAlleles());
            }
        }
        return alleles;
    }

    public void print(PrintWriter str) {
        for (int c=0; c<popdata.chromCount(); c++) {
            for (int h=0; h<getPloidy(); h++) {
                str.println(c+" "+h+" : "+haplostruct[c][h]);
            }
        }
    }
    
    /**
     * recPositions lists for a given chromosome all recombination positions
     * over all haplotypes for that chromosome.
     * Identical positions are counted only once.
     * For gametes this produces all recombination positions
     * 
     * @return ordered array of all different recombination positions
     */
    double[] recPositions(int chrom) {
        //get all recPositions of the first haplotype:
        double[] recpos = new double[0];
        for (int p=0; p<getPloidy(); p++) {
            ArrayList<Double> hrecs = haplostruct[chrom][p].getRecombPos();
            int rp=0; //index to recpos
            for (Double hrec : hrecs) {
                //first: skip all smaller recpos items:
                while (rp<recpos.length && recpos[rp] < hrec) {
                    rp++;
                }
                //if hrecs value not in recpos: insert
                if (rp>=recpos.length || recpos[rp] != hrec) {
                    double[] tmp = new double[recpos.length+1];
                    System.arraycopy(recpos, 0, tmp, 0, rp);
                    tmp[rp] = hrec;
                    System.arraycopy(recpos, rp, tmp, rp + 1, recpos.length - rp);
                    recpos = tmp;
                }
                //else hrec already in recpos, do nothing
                //rc (still) points to recpos value identical to hrecs value
            }    
        }
        assert recpos[0]==popdata.getChrom(chrom).getHeadPos();
        //delete the first (headPos) position:
        double[] result = new double[recpos.length-1];
        for (int i=0; i<result.length; i++) {
            result[i] = recpos[i+1];
        }
        return result;
    }

}
