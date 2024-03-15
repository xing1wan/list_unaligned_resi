from pymol import cmd, stored

def list_unaligned_resi(ref_structure, tgt_structure, cutoff=2.0):
    """
    Author: Xing Wan Date: 2024-03-15
    Inspired by rmsdByRes Zhenting Gao on 7/28/2016 and rmsdCA Yufeng Tong 2020-01-29
    
    Workflow
     Load reference and target protein stuctures
     1. Select alpha carbons of the two structures
     2. Align those CAs
     3. Calculate the atomatic distances
     4. Select aligned and unaligned residues
    """
    
    stored.ca_ref_crd = []
    stored.ca_target_crd = []
    stored.all_ref_residues = set()
    stored.all_tgt_residues = set()
    stored.aligned_residues_ref = set()
    stored.aligned_residues_tgt = set()

    # Generate dynamic selection names for alpha carbons in reference and target structures
    ref_ca = "{}_ca".format(ref_structure)
    tgt_ca = "{}_ca".format(tgt_structure)

    # Create selections for alpha carbons
    cmd.select(ref_ca, "({}) and name CA".format(ref_structure))
    cmd.select(tgt_ca, "({}) and name CA".format(tgt_structure))
    
    # Perform alignment based on alpha carbon selections, set cutoff the same as global cutoff
    alignment = cmd.align(tgt_ca, ref_ca, cutoff = 2)
    print("Alignment result:", alignment)
    
    # Get coordinates of alpha carbons post-alignment
    cmd.iterate_state(1, ref_ca, "stored.ca_ref_crd.append((resi, resn, (x,y,z))); stored.all_ref_residues.add(resi)")
    cmd.iterate_state(1, tgt_ca, "stored.ca_target_crd.append((resi, resn, (x,y,z))); stored.all_tgt_residues.add(resi)")
    
    def distance(coord1, coord2):
        """Calculate the Euclidean distance between two points."""
        return ((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)**0.5
    
    aligned_pairs = []
    for ref in stored.ca_ref_crd:
        closest_distance = None
        closest_target = None
        for tgt in stored.ca_target_crd:
            dist = distance(ref[2], tgt[2])
            if closest_distance is None or dist < closest_distance:
                closest_distance = dist
                closest_target = tgt
                
        if closest_distance < cutoff:
            print("Aligned: Ref {}{} - Target {}{} (Distance: {:.2f} Å)".format(ref[1], ref[0], closest_target[1], closest_target[0], closest_distance))
            aligned_pairs.append((ref[0], closest_target[0]))
            stored.aligned_residues_ref.add(ref[0])
            stored.aligned_residues_tgt.add(closest_target[0])
    
    # Identifying unaligned residues
    unaligned_ref = stored.all_ref_residues - stored.aligned_residues_ref
    unaligned_tgt = stored.all_tgt_residues - stored.aligned_residues_tgt
    aligned_ref = stored.aligned_residues_ref
    aligned_tgt = stored.aligned_residues_tgt
    
    # Creating selections for unaligned residues
    if unaligned_ref:
        cmd.select("unaligned_ref", " or ".join(["{} and resi {}".format(ref_structure, resi) for resi in unaligned_ref]))
        #cmd.show("surface", "unaligned_ref")
        cmd.color("red", "unaligned_ref")
    if unaligned_tgt:
        cmd.select("unaligned_tgt", " or ".join(["{} and resi {}".format(tgt_structure, resi) for resi in unaligned_tgt]))
        # cmd.show("surface", "unaligned_tgt")
        cmd.color("blue", "unaligned_tgt")
    # Creating selections for aligned residues
    if aligned_ref:
        cmd.select("aligned_ref", " or ".join(["{} and resi {}".format(ref_structure, resi) for resi in aligned_ref]))
        # cmd.show("surface", "aligned_ref")
        cmd.color("grey80", "aligned_ref")
    if aligned_tgt:
        cmd.select("aligned_tgt", " or ".join(["{} and resi {}".format(tgt_structure, resi) for resi in aligned_tgt]))
        # cmd.show("surface", "aligned_tgt")
        cmd.color("grey80", "aligned_tgt")


    # print("Unaligned Reference Residues:", unaligned_ref)
    # print("Unaligned Target Residues:", unaligned_tgt)

# Example usage:
# list_unaligned_resi('/ref_structure//A', '/tgt_structure//A', 1.0)
print("""Example usage:
list_unaligned_resi('/ref_structure//A', '/tgt_structure//A', 1.0)""")
