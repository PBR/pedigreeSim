/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

/**
 * Bivalent extends the abstract class Multivalent.
 * Actually the only addition is a constructor that checks if the 
 * supplied haplostruct is valid for a bivalent.
 * A Bivalent is composed of two homologous chromosomes (4 chromatids).
 * The order of the chromosomes is randomized when constructing the Bivalent.
 * @author Roeland Voorrips
 */
public class Bivalent extends Multivalent {

    /**
     * The order of the chromosomes is NOT randomized in the constructor!
     * Instead it is determined by the caller (Individual.doMeiosis)
     * (because for a Quadrivalent the randomization must be
     * done there)
     * @param haplostruct
     * @throws Exception
     */
    public Bivalent(HaploStruct[] haplostruct)
            throws Exception {
        super();
        if (haplostruct==null || haplostruct.length!=2 ||
                haplostruct[0]==null || haplostruct[1]==null ||
                !haplostruct[0].getChrom().equals(haplostruct[1].getChrom())) {
            throw new Exception("Error in Bivalent constructor parameters");
        }
        this.chrom = haplostruct[0].getChrom();
        this.popdata = chrom.getPopdata();
        this.tools = popdata.tools;
        this.rand = tools.rand;
        this.haplostruct = haplostruct;
    }

    public Bivalent(HaploStruct haplostruct0, HaploStruct haplostruct1)
            throws Exception {
        this(new HaploStruct[] {haplostruct0,haplostruct1} );
    }


}
