#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Tianyu (Sky) Lu (tianyu@lu.fm)

from core import SequenceTree


def test_multiple_cuts():
    """
    Test reconstruction algorithm for multiple cutting scenarios
    Create a scenario where multiple donors cut a sequence, then perform reconstruction 
    and verify that all cutting relationships are properly processed.
    Returns:
        tuple: (List of reconstructed donor records, success status)
    """
    # Initialize tree data structure
    tree = SequenceTree("ATGCATGCATGCATGCATGCATGCATGCATGCATGC")  # Original sequence

    # Create a simple scenario with multiple cuts
    # The original sequence is cut by three different donors:
    # - donor 1 cuts at position 10
    # - donor 2 cuts at position 20
    # - donor 3 cuts at position 30
    # Using real DNA sequences
    first_cutter = "GTACGTAC"
    second_cutter = "CCGGAATT"
    third_cutter = "TTAGGCCA"
    
    # Insert donor sequences at respective positions
    tree.insert(10, first_cutter, "cutter1")
    tree.insert(20, second_cutter, "cutter2")
    tree.insert(30, third_cutter, "cutter3")
    
    # Reconstruct donors
    _, reconstructed_records = tree.donors("test")
    
    # Define expected results - explicitly based on expected behavior of SequenceTree
    expected_results = {
        "clean_recon_count": 3,  # Expect 3 clean reconstructions (one for each cut point)
        "full_recon_count": 3,   # Expect 3 full reconstructions (one for each cutter)
        "expected_cutters": [first_cutter, second_cutter, third_cutter],
        "min_cutter_matches": 2  # At least 2 cutters should be found in full reconstructions
    }
    
    # Verify correctness
    success = True
    error_messages = []
    
    # Check clean reconstruction count
    clean_recon = [rec for rec in reconstructed_records if rec.annotations.get("reconstruction_type") == "clean"]
    if len(clean_recon) != expected_results["clean_recon_count"]:
        error_messages.append(f"Error: Should have {expected_results['clean_recon_count']} clean reconstructions, but actually have {len(clean_recon)}")
        success = False
    
    # Check full reconstruction count
    full_recon = [rec for rec in reconstructed_records if rec.annotations.get("reconstruction_type") == "full"]
    if len(full_recon) != expected_results["full_recon_count"]:
        error_messages.append(f"Error: Should have {expected_results['full_recon_count']} full reconstructions, but actually have {len(full_recon)}")
        success = False
    
    # Validate full reconstruction sequence content - check if it contains cutter sequences
    found_cutters = []
    for cutter in expected_results["expected_cutters"]:
        found = False
        for rec in full_recon:
            seq_str = str(rec.seq)
            # Check complete or partial sequence
            if cutter in seq_str or cutter[:4] in seq_str:
                found = True
                found_cutters.append(cutter)
                break
    
    if len(found_cutters) < expected_results["min_cutter_matches"]:
        error_messages.append(f"Error: Should find at least {expected_results['min_cutter_matches']} cutter sequences, but only found {len(found_cutters)}")
        success = False
    
    # Print test results
    if success:
        print("Multiple cutting test successful! All cutting relationships were correctly processed.")
    else:
        print("Multiple cutting test failed! Please check the reconstruction algorithm.")
        for error in error_messages:
            print(f"  {error}")
    
    return reconstructed_records, success


def test_comprehensive_nesting():
    """
    Comprehensive test of various complex nesting and cutting scenarios to verify 
    the correctness of the reconstruction algorithm.
    Tests include multiple scenarios:
    - Single insertion
    - Adjacent position insertions
    - Nested insertions
    - Multiple nesting
    - Cutting donors
    - Multiple cuts
    - Chain cutting
    - Boundary cutting
    - Random multi-donor complex networks
    - Empty sequence handling
    - Extreme length sequences
    - Completely overlapping cuts
    Returns:
        bool: Whether all tests passed
    """
    test_results = []
    
    # ======== Scenario 1: Single donor insertion ========
    print("\n===== Testing Scenario 1: Single donor insertion =====")
    # Create initial sequence
    tree = SequenceTree("ATGC")
    # Insert donor
    tree.insert(2, "TTT", "donor1")
    
    # Define expected results
    expected = {
        "total_reconstructions": 2,         # Expected total of 2 reconstructions
        "clean_reconstructions": 1,         # Expected 1 clean reconstruction
        "full_reconstructions": 1,          # Expected 1 full reconstruction
        "donor_sequence_present": "TTT"     # Expected this sequence to be in reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    # Check total reconstruction count
    if reconstructed is None or len(reconstructed) != expected["total_reconstructions"]:
        success = False
        reason = f"Expected {expected['total_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    else:
        # Check clean reconstruction count
        clean_count = len([r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"])
        if clean_count != expected["clean_reconstructions"]:
            success = False
            reason = f"Expected {expected['clean_reconstructions']} clean reconstructions, but got {clean_count}"
        
        # Check full reconstruction count
        full_count = len([r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"])
        if full_count != expected["full_reconstructions"]:
            success = False
            reason = f"Expected {expected['full_reconstructions']} full reconstructions, but got {full_count}"
        
        # Check sequence content
        if not any(expected["donor_sequence_present"] in str(r.seq) for r in reconstructed):
            success = False
            reason = f"Expected sequence {expected['donor_sequence_present']} not found in reconstructions"
    
    # Record test results
    test_results.append(("Scenario 1", success, reason if not success else "Single donor insertion, created expected reconstructions"))
    if success:
        print("✓ Scenario 1 test passed: Single donor insertion, created expected reconstructions")
    else:
        print(f"✗ Scenario 1 test failed: {reason}")
    
    # ======== Scenario 2: Adjacent double insertion ========
    print("\n===== Testing Scenario 2: Adjacent double insertion =====")
    # Create initial sequence
    tree = SequenceTree("ATGC")
    # Insert donor1
    tree.insert(2, "TTT", "donor1")
    # Insert donor2 (adjacent position)
    tree.insert(5, "GGG", "donor2")
    
    # Define expected results
    expected = {
        "total_reconstructions": 2,  # Expected total of 2 reconstructions
        "donor_sequences": ["TTT", "GGG"]  # Expected to contain these two sequences
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) != expected["total_reconstructions"]:
        success = False
        reason = f"Expected {expected['total_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    else:
        # Check sequence content
        for donor_seq in expected["donor_sequences"]:
            if not any(donor_seq in str(r.seq) for r in reconstructed):
                success = False
                reason = f"Expected sequence {donor_seq} not found in reconstructions"
                break
    
    # Record test results
    test_results.append(("Scenario 2", success, reason if not success else "Adjacent donors, created expected reconstructions"))
    if success:
        print("✓ Scenario 2 test passed: Adjacent donors, created expected reconstructions")
    else:
        print(f"✗ Scenario 2 test failed: {reason}")
    
    # ======== Scenario 3: Nested insertion ========
    print("\n===== Testing Scenario 3: Nested insertion =====")
    # Create initial sequence
    tree = SequenceTree("ATGC")
    # Insert donor1
    tree.insert(2, "TTTAAA", "donor1")
    # Insert donor2(in donor1)
    tree.insert(5, "GGG", "donor2")
    
    # Define expected results
    expected = {
        "min_reconstructions": 2,  # Expected at least 2 reconstructions
        "min_full_reconstructions": 1,  # Expected at least 1 full reconstruction
        "should_contain_nested": True,  # Expected at least one reconstruction to contain nested sequence
        "outer_donor": "TTTAAA",
        "inner_donor": "GGG"
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    else:
        # Check full reconstruction count
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        if len(full_recon) < expected["min_full_reconstructions"]:
            success = False
            reason = f"Expected at least {expected['min_full_reconstructions']} full reconstructions, but got {len(full_recon)}"
        else:
            # Check nested sequence
            found_nested = False
            for rec in full_recon:
                seq_str = str(rec.seq)
                if expected["inner_donor"] in seq_str and (expected["outer_donor"] in seq_str or 
                       expected["outer_donor"][:3] in seq_str or expected["outer_donor"][3:] in seq_str):
                    found_nested = True
                    break
            
            if expected["should_contain_nested"] and not found_nested:
                success = False
                reason = "No reconstruction containing expected nested sequence found"
    
    # Record test results
    test_results.append(("Scenario 3", success, reason if not success else "Nested insertion reconstruction correct"))
    if success:
        print("✓ Scenario 3 test passed: Nested insertion reconstruction correct")
    else:
        print(f"✗ Scenario 3 test failed: {reason}")
        
    # ======== Scenario 4: Multiple nesting insertion ========
    print("\n===== Testing Scenario 4: Multiple nesting insertion =====")
    # Create initial sequence
    tree = SequenceTree("ATGC")
    # Insert donor1
    tree.insert(2, "TTTAAA", "donor1")
    # Insert donor2(in donor1)
    tree.insert(5, "GGCCC", "donor2")
    # Insert donor3(in donor2)
    tree.insert(7, "AAA", "donor3")
    
    # Define expected results
    expected = {
        "min_reconstructions": 2,  # Expected at least 2 reconstructions
        "min_nested_sequences": 2,  # Expected at least 2 layers of nesting
        "donor_sequences": ["TTTAAA", "GGCCC", "AAA"]  # Possible donor sequences
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    else:
        # Get all full reconstructions
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        
        # Check if at least one reconstruction contains multiple layers of nesting
        max_nested_count = 0
        for rec in full_recon:
            seq_str = str(rec.seq)
            nested_count = 0
            for donor in expected["donor_sequences"]:
                if donor in seq_str or (len(donor) >= 6 and (donor[:3] in seq_str or donor[3:] in seq_str)):
                    nested_count += 1
            max_nested_count = max(max_nested_count, nested_count)
        
        if max_nested_count < expected["min_nested_sequences"]:
            success = False
            reason = f"Expected to find at least {expected['min_nested_sequences']} layers of nesting, but found {max_nested_count} layers"
    
    # Record test results
    test_results.append(("Scenario 4", success, reason if not success else "Multiple nesting insertion reconstruction correct"))
    if success:
        print("✓ Scenario 4 test passed: Multiple nesting insertion reconstruction correct")
    else:
        print(f"✗ Scenario 4 test failed: {reason}")
    
    # ======== Scenario 5: Cutting donor ========
    print("\n===== Testing Scenario 5: Cutting donor ========")
    # Create initial sequence
    tree = SequenceTree("ATGC")
    # Insert donor1
    tree.insert(2, "TTAAA", "donor1")
    # Insert donor2(cut off donor1)
    tree.insert(4, "GGG", "donor2")
    
    # Define expected results
    expected = {
        "min_reconstructions": 2,  # Expected at least 2 reconstructions
        "has_clean_reconstructions": True,  # Expected clean reconstruction
        "has_full_reconstructions": True,  # Expected full reconstruction
        "cutter_sequence": "GGG"  # Expected to find cutter sequence in full reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    else:
        # Verify full and clean reconstruction
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        
        if expected["has_clean_reconstructions"] and not clean_recon:
            success = False
            reason = "Missing clean reconstruction"
        elif expected["has_full_reconstructions"] and not full_recon:
            success = False
            reason = "Missing full reconstruction"
        else:
            # Check if at least one full reconstruction contains donor2
            has_donor2 = any(expected["cutter_sequence"] in str(r.seq) for r in full_recon)
            
            if not has_donor2:
                success = False
                reason = f"No full reconstruction containing {expected['cutter_sequence']}"
    
    # Record test results
    test_results.append(("Scenario 5", success, reason if not success else "Cutting donor reconstruction correct"))
    if success:
        print("✓ Scenario 5 test passed: Cutting donor reconstruction correct")
    else:
        print(f"✗ Scenario 5 test failed: {reason}")
    
    # ======== Scenario 6: Multiple cuts ========
    print("\n===== Testing Scenario 6: Multiple cuts =====")
    # Create initial sequence
    tree = SequenceTree("ATGC")
    # Insert donor1
    tree.insert(2, "TTTAAACCC", "donor1")
    # Insert donor2(cut off donor1)
    tree.insert(4, "GGG", "donor2")
    # Insert donor3(cut off donor1 again)
    tree.insert(8, "TTT", "donor3")
    
    # Define expected results
    expected = {
        "min_reconstructions": 1,  # Expected at least 1 reconstruction
        "has_clean_reconstructions": True,  # Expected clean reconstruction
        "has_full_reconstructions": True,  # Expected full reconstruction
        "cutter_sequences": ["GGG", "TTT"]  # Expected to find at least one cutter sequence in full reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    else:
        # Verify clean and full reconstruction existence
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        
        if expected["has_clean_reconstructions"] and not clean_recon:
            success = False
            reason = "Missing clean reconstruction"
        elif expected["has_full_reconstructions"] and not full_recon:
            success = False
            reason = "Missing full reconstruction"
        else:
            # Check if full reconstruction contains any cutter sequence
            has_cutter = False
            for cutter in expected["cutter_sequences"]:
                if any(cutter in str(r.seq) for r in full_recon):
                    has_cutter = True
                    break
            
            if not has_cutter:
                success = False
                reason = "No full reconstruction containing cutter sequence"
    
    # Record test results
    test_results.append(("Scenario 6", success, reason if not success else "Multiple cuts reconstruction correct"))
    if success:
        print("✓ Scenario 6 test passed: Multiple cuts reconstruction correct")
    else:
        print(f"✗ Scenario 6 test failed: {reason}")
    
    # ======== Scenario 7: Chain cutting ========
    print("\n===== Testing Scenario 7: Chain cutting =====")
    # Create initial sequence
    tree = SequenceTree("ATGC")
    # Insert donor1
    tree.insert(2, "TTTAAA", "donor1")
    # Insert donor2(cut off donor1)
    tree.insert(4, "GGGCCC", "donor2")
    # Insert donor3(cut off donor2)
    tree.insert(6, "AAA", "donor3")
    
    # Define expected results
    expected = {
        "min_reconstructions": 3,  # Expected at least 3 reconstructions
        "min_full_reconstructions": 1,  # Expected at least 1 full reconstruction
        "min_clean_reconstructions": 1,  # Expected at least 1 clean reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    else:
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        
        if len(full_recon) < expected["min_full_reconstructions"]:
            success = False
            reason = f"Expected at least {expected['min_full_reconstructions']} full reconstructions, but got {len(full_recon)}"
        elif len(clean_recon) < expected["min_clean_reconstructions"]:
            success = False
            reason = f"Expected at least {expected['min_clean_reconstructions']} clean reconstructions, but got {len(clean_recon)}"
    
    # Record test results
    test_results.append(("Scenario 7", success, reason if not success else "Chain cutting reconstruction correct"))
    if success:
        print("✓ Scenario 7 test passed: Chain cutting reconstruction correct")
    else:
        print(f"✗ Scenario 7 test failed: {reason}")
    
    # ======== Scenario 8: First and last insertion boundary cases ========
    print("\n===== Testing Scenario 8: First and last insertion boundary cases =====")
    # Create initial sequence
    tree = SequenceTree("ATGC")
    # Insert donor at start of sequence
    tree.insert(1, "TTT", "donor1")
    # Insert donor at end of sequence
    tree.insert(7, "GGG", "donor2")
    
    # Define expected results
    expected = {
        "min_reconstructions": 1  # Expected at least 1 reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    
    # Record test results
    test_results.append(("Scenario 8", success, reason if not success else "First and last insertion boundary cases correct processing"))
    if success:
        print("✓ Scenario 8 test passed: First and last insertion boundary cases correct processing")
    else:
        print(f"✗ Scenario 8 test failed: {reason}")
    
    # ======== Scenario 9: Multiple insertions at same position ========
    print("\n===== Testing Scenario 9: Multiple insertions at same position =====")
    # Create initial sequence
    tree = SequenceTree("ATGC")
    # Insert two donors at same position
    tree.insert(2, "TTT", "donor1")
    tree.insert(2, "GGG", "donor2")
    
    # Define expected results
    expected = {
        "min_reconstructions": 1  # Expected at least 1 reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    
    # Record test results
    test_results.append(("Scenario 9", success, reason if not success else "Multiple insertions at same position correct processing"))
    if success:
        print("✓ Scenario 9 test passed: Multiple insertions at same position correct processing")
    else:
        print(f"✗ Scenario 9 test failed: {reason}")
    
    # ======== Scenario 10: Random multi-donor complex network ========
    print("\n===== Testing Scenario 10: Random multi-donor complex network =====")
    # Create initial sequence
    tree = SequenceTree("ATGCATGCATGC")
    # Insert multiple mutually cutting donors
    tree.insert(3, "AAATTT", "donor1")
    tree.insert(5, "CCCGGG", "donor2")
    tree.insert(8, "TTTAAA", "donor3")
    tree.insert(10, "GGGCCC", "donor4")
    
    # Define expected results
    expected = {
        "min_reconstructions": 4  # Expected at least 4 reconstructions
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    
    # Record test results
    test_results.append(("Scenario 10", success, reason if not success else "Random multi-donor complex network correct processing"))
    if success:
        print("✓ Scenario 10 test passed: Random multi-donor complex network correct processing")
    else:
        print(f"✗ Scenario 10 test failed: {reason}")
    
    # ======== Scenario 11: Empty sequence handling ========
    print("\n===== Testing Scenario 11: Empty sequence handling =====")
    # Create empty sequence
    tree = SequenceTree("")
    # Insert donor to empty sequence
    tree.insert(1, "TTT", "donor1")
    
    # Define expected results - based on actual behavior of SequenceTree
    expected = {
        "should_have_reconstructions": False  # Expected no reconstructions (because no nesting or cutting)
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if expected["should_have_reconstructions"] == False:
        if reconstructed and len(reconstructed) > 0:
            success = False
            reason = f"Expected no reconstructions, but got {len(reconstructed)} reconstructions"
    else:
        if reconstructed is None or len(reconstructed) == 0:
            success = False
            reason = "Expected reconstructions, but none obtained"
    
    # Record test results
    test_results.append(("Scenario 11", success, reason if not success else "Empty sequence correct processing"))
    if success:
        print("✓ Scenario 11 test passed: Empty sequence correct processing")
    else:
        print(f"✗ Scenario 11 test failed: {reason}")
    
    # ======== Scenario 12: Extremely long sequence handling ========
    print("\n===== Testing Scenario 12: Extremely long sequence handling =====")
    # Create a longer sequence (1000 bases)
    long_seq = "A" * 400 + "T" * 300 + "G" * 200 + "C" * 100
    tree = SequenceTree(long_seq)
    # Insert two mutually nested donors
    tree.insert(500, "AAATTT", "donor1")
    tree.insert(502, "GGG", "donor2")
    
    # Define expected results
    expected = {
        "min_reconstructions": 1,  # Expected at least 1 reconstruction
        "min_full_reconstructions": 1  # Expected at least 1 full reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    else:
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        if len(full_recon) < expected["min_full_reconstructions"]:
            success = False
            reason = f"Expected at least {expected['min_full_reconstructions']} full reconstructions, but got {len(full_recon)}"
    
    # Record test results
    test_results.append(("Scenario 12", success, reason if not success else "Extremely long sequence correct processing"))
    if success:
        print("✓ Scenario 12 test passed: Extremely long sequence correct processing")
    else:
        print(f"✗ Scenario 12 test failed: {reason}")
    
    # ======== Scenario 13: Completely overlapping cuts ========
    print("\n===== Testing Scenario 13: Completely overlapping cuts =====")
    # Create initial sequence
    tree = SequenceTree("ATGCATGC")
    # Create a donor
    tree.insert(2, "AAATTTCCC", "donor1")
    # Create two donors simultaneously cutting the same position
    tree.insert(5, "GGG", "cutter1")
    tree.insert(5, "TTT", "cutter2")
    
    # Define expected results
    expected = {
        "min_reconstructions": 1,  # Expected at least 1 reconstruction
        "min_clean_reconstructions": 1,  # Expected at least 1 clean reconstruction
        "min_full_reconstructions": 1   # Expected at least 1 full reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    else:
        # Should have at least one clean reconstruction and one full reconstruction
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        
        if len(clean_recon) < expected["min_clean_reconstructions"]:
            success = False
            reason = f"Expected at least {expected['min_clean_reconstructions']} clean reconstructions, but only got {len(clean_recon)}"
        elif len(full_recon) < expected["min_full_reconstructions"]:
            success = False
            reason = f"Expected at least {expected['min_full_reconstructions']} full reconstructions, but only got {len(full_recon)}"
    
    # Record test results
    test_results.append(("Scenario 13", success, reason if not success else "Completely overlapping cuts correct processing"))
    if success:
        print("✓ Scenario 13 test passed: Completely overlapping cuts correct processing")
    else:
        print(f"✗ Scenario 13 test failed: {reason}")
    
    # ======== Scenario 14: Boundary cutting ========
    print("\n===== Testing Scenario 14: Boundary cutting =====")
    # Create initial sequence
    tree = SequenceTree("ATGCATGC")
    # Create a donor
    tree.insert(2, "AAATTTCCC", "donor1")
    # Create a cutter cutting at donor boundary
    tree.insert(2, "GGG", "cutter1")  # Cutting at donor start
    
    # Define expected results
    expected = {
        "min_reconstructions": 1  # Expected at least 1 reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    
    # Record test results
    test_results.append(("Scenario 14", success, reason if not success else "Boundary cutting correct processing"))
    if success:
        print("✓ Scenario 14 test passed: Boundary cutting correct processing")
    else:
        print(f"✗ Scenario 14 test failed: {reason}")
    
    # ======== Scenario 15: Multiple donors simultaneously cutting another donor at different positions ========
    print("\n===== Testing Scenario 15: Multiple donors simultaneously cutting another donor at different positions =====")
    # Create initial sequence
    tree = SequenceTree("ATGCATGC")
    # Create main donor
    tree.insert(2, "AAAAAAAAAATTTTTTTTTTGGGGGGGGGG", "main_donor")  # 30bp
    # Create three cutters
    tree.insert(5, "CCC", "cutter1")  # Cutting front segment
    tree.insert(15, "TTT", "cutter2")  # Cutting middle segment
    tree.insert(25, "GGG", "cutter3")  # Cutting back segment
    
    # Define expected results
    expected = {
        "min_reconstructions": 1,  # Expected at least 1 reconstruction
        "has_clean_reconstructions": True  # Expected clean reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    else:
        # Check clean reconstruction
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        if expected["has_clean_reconstructions"] and not clean_recon:
            success = False
            reason = "Missing main donor clean reconstruction"
    
    # Record test results
    test_results.append(("Scenario 15", success, reason if not success else "Multiple donors simultaneously cutting another donor at different positions correct processing"))
    if success:
        print("✓ Scenario 15 test passed: Multiple donors simultaneously cutting another donor at different positions correct processing")
        # Check full reconstruction count
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        print(f"   Full reconstruction count: {len(full_recon)}")
    else:
        print(f"✗ Scenario 15 test failed: {reason}")
    
    # ======== Scenario 16: Complex cutting network - Multiple donors mutually cutting each other ========
    print("\n===== Testing Scenario 16: Complex cutting network - Multiple donors mutually cutting each other =====")
    # Create initial sequence
    tree = SequenceTree("ATGCATGCATGCATGCATGC")
    # Create 5 donors, forming complex mutually cutting network
    tree.insert(2, "AAAAAAAAAA", "donor_a")   # A
    tree.insert(12, "TTTTTTTTTT", "donor_b")  # B
    tree.insert(22, "GGGGGGGGGG", "donor_c")  # C
    tree.insert(32, "CCCCCCCCCC", "donor_d")  # D
    tree.insert(42, "ATATATATAT", "donor_e")  # E
    
    # Define expected results
    expected = {
        "min_reconstructions": 1  # Expected at least 1 reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    
    # Record test results
    test_results.append(("Scenario 16", success, reason if not success else f"Complex cutting network correct processing, generated {len(reconstructed) if reconstructed else 0} reconstructions"))
    if success:
        print(f"✓ Scenario 16 test passed: Complex cutting network correct processing, generated {len(reconstructed) if reconstructed else 0} reconstructions")
    else:
        print(f"✗ Scenario 16 test failed: {reason}")
    
    # ======== Scenario 17: Sequence end boundary cutting and first and last tail special cutting ========
    print("\n===== Testing Scenario 17: Sequence end boundary cutting and first and last tail special cutting =====")
    # Create initial sequence
    tree = SequenceTree("ATGCATGC")
    # Create donor sequence
    tree.insert(2, "AAAATTTTGGGG", "donor")  # 12bp
    # Create donor cutting at first and last positions
    tree.insert(2, "CCC", "cutter1")  # Cutting at donor start
    tree.insert(14, "TTT", "cutter2")  # Cutting at donor end
    
    # Define expected results
    expected = {
        "min_reconstructions": 1,  # Expected at least 1 reconstruction
        "has_clean_reconstructions": True  # Expected clean reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    elif expected["has_clean_reconstructions"]:
        # Check if correctly processed multiple cuts
        clean_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "clean"]
        if not clean_recon:
            success = False
            reason = "Missing donor clean reconstruction"
    
    # Record test results
    test_results.append(("Scenario 17", success, reason if not success else "Sequence end boundary cutting and first and last tail special cutting correct processing"))
    if success:
        print("✓ Scenario 17 test passed: Sequence end boundary cutting and first and last tail special cutting correct processing")
    else:
        print(f"✗ Scenario 17 test failed: {reason}")
    
    # ======== Scenario 18: Deep nesting ========
    print("\n===== Testing Scenario 18: Deep nesting =====")
    # Create initial sequence
    tree = SequenceTree("ATGCATGC")
    # Create 5 layers of nested donors
    tree.insert(2, "AAAAAAAA", "donor1")  # Outer layer
    tree.insert(4, "TTTTTTTT", "donor2")  # 2nd layer
    tree.insert(6, "GGGGGGGG", "donor3")  # 3rd layer
    tree.insert(8, "CCCCCCCC", "donor4")  # 4th layer
    tree.insert(10, "AAAATTTT", "donor5") # Inner layer
    
    # Define expected results
    expected = {
        "min_reconstructions": 1,  # Expected at least 1 reconstruction
        "min_full_reconstructions": 1  # Expected at least 1 full reconstruction
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if reconstructed is None or len(reconstructed) < expected["min_reconstructions"]:
        success = False
        reason = f"Expected at least {expected['min_reconstructions']} reconstructions, but got {len(reconstructed) if reconstructed else 0} reconstructions"
    else:
        full_recon = [r for r in reconstructed if r.annotations.get("reconstruction_type") == "full"]
        if len(full_recon) < expected["min_full_reconstructions"]:
            success = False
            reason = f"Expected at least {expected['min_full_reconstructions']} full reconstructions, but got {len(full_recon)}"
    
    # Record test results
    test_results.append(("Scenario 18", success, reason if not success else "Deep nesting correct processing"))
    if success:
        print("✓ Scenario 18 test passed: Deep nesting correct processing")
        print(f"   Generated {len(full_recon)} full reconstructions")
    else:
        print(f"✗ Scenario 18 test failed: {reason}")
    
    # ======== Scenario 19: Extremely short sequence and empty sequence cutting processing ========
    print("\n===== Testing Scenario 19: Extremely short sequence and empty sequence cutting processing =====")
    # Create initial sequence
    tree = SequenceTree("ATGC")
    # Create extremely short donor sequence
    tree.insert(2, "A", "short_donor")
    # Create a donor cutting extremely short sequence
    tree.insert(2, "TTT", "cutter")
    
    # Define expected results
    expected = {
        "should_have_reconstructions": True  # Expected reconstructions
    }
    
    # Verify results
    _, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    if expected["should_have_reconstructions"]:
        if reconstructed is None or len(reconstructed) == 0:
            success = False
            reason = "Expected reconstructions, but none obtained"
    else:
        if reconstructed and len(reconstructed) > 0:
            success = False
            reason = f"Expected no reconstructions, but got {len(reconstructed)} reconstructions"
    
    # Record test results
    test_results.append(("Scenario 19", success, reason if not success else "Extremely short sequence correct processing"))
    if success:
        print("✓ Scenario 19 test passed: Extremely short sequence correct processing")
    else:
        print(f"✗ Scenario 19 test failed: {reason}")
    
    # ======== Display final test results ========
    print("\n===== Comprehensive test results =====")
    passed = sum(1 for _, result, _ in test_results if result)
    total = len(test_results)
    print(f"Test completed: {passed}/{total} passed")
    if passed == total:
        print("✓ All scenarios passed! Reconstruction algorithm works correctly")
    else:
        print("✗ Some scenarios failed, please check the reconstruction algorithm")
        for name, result, message in test_results:
            if not result:
                print(f"  - {name}: {message}")
    
    return passed == total


def test_multiple_cuts_fragments_distinction():
    """
    Test improvement of fragment distinction function for multiple cutting scenarios
    Verify if it can correctly distinguish different fragments produced by different cutters
    
    Returns:
        bool: Whether test passed
    """
    print("\n===== Testing multiple cutting fragment distinction function =====")
    
    # Initialize tree data structure
    tree = SequenceTree("ATGCATGCATGCATGCATGCATGCATGCATGCATGC")  # Original sequence

    # Create a scenario with multiple cuts
    # The original sequence is cut by three different donors:
    # - donor 1 cuts at position 10
    # - donor 2 cuts at position 20
    # - donor 3 cuts at position 30
    
    # Insert donor sequences at respective positions
    tree.insert(10, "GTACGTAC", "cutter1")
    tree.insert(20, "CCGGAATT", "cutter2")
    tree.insert(30, "TTAGGCCA", "cutter3")
    
    # Get event journal
    event_journal = tree.event_journal
    
    # Define expected results
    expected = {
        "event_count": 3,  # Expected 3 event records
        "donor_uids": [],  # Will be filled during test
        "fragment_uids": [],  # Will be filled during test
        "has_valid_fragment_info": True,  # Expected all fragments to have valid information
        "fragment_cutter_match": True,  # Expected cutter information in fragment to match donor_uid in event
    }
    
    # Test fragment distinction function
    success = True
    error_messages = []
    
    # 1. Verify if event records are complete
    if len(event_journal.events) != expected["event_count"]:
        error_messages.append(f"Error: Should record {expected['event_count']} events, but found {len(event_journal.events)}")
        success = False
    
    # 2. Get event information and fill expected results UID list
    for event in event_journal.events:
        print(f"Event {event.event_id}: donor({event.donor_uid}) → target({event.target_uid}) → [L({event.left_uid}), R({event.right_uid})]")
        
        # Record donor and fragment UID
        expected["donor_uids"].append(event.donor_uid)
        expected["fragment_uids"].append(event.left_uid)
        expected["fragment_uids"].append(event.right_uid)
        
        # Check donor is correctly marked
        if not event_journal.is_donor(event.donor_uid):
            error_messages.append(f"Error: UID {event.donor_uid} should be marked as donor")
            success = False
        
        # Check fragment is correctly marked
        if not event_journal.is_fragment(event.left_uid):
            error_messages.append(f"Error: UID {event.left_uid} should be marked as fragment")
            success = False
        
        if not event_journal.is_fragment(event.right_uid):
            error_messages.append(f"Error: UID {event.right_uid} should be marked as fragment")
            success = False
        
        # Check fragment information
        left_info = event_journal.get_fragment_info(event.left_uid)
        right_info = event_journal.get_fragment_info(event.right_uid)
        
        if left_info:
            orig_uid, is_left, cutter_uid = left_info
            print(f"   Left fragment {event.left_uid}: From {orig_uid}, "
                  f"{'Left' if is_left else 'Right'} side, Cutter {cutter_uid}")
            
            # Verify cutter information
            if expected["fragment_cutter_match"] and cutter_uid != event.donor_uid:
                error_messages.append(f"Error: Left fragment {event.left_uid} cutter should be {event.donor_uid}, but got {cutter_uid}")
                success = False
        elif expected["has_valid_fragment_info"]:
            error_messages.append(f"Error: Left fragment {event.left_uid} information not found")
            success = False
        
        if right_info:
            orig_uid, is_left, cutter_uid = right_info
            print(f"   Right fragment {event.right_uid}: From {orig_uid}, "
                  f"{'Left' if is_left else 'Right'} side, Cutter {cutter_uid}")
            
            # Verify cutter information
            if expected["fragment_cutter_match"] and cutter_uid != event.donor_uid:
                error_messages.append(f"Error: Right fragment {event.right_uid} cutter should be {event.donor_uid}, but got {cutter_uid}")
                success = False
        elif expected["has_valid_fragment_info"]:
            error_messages.append(f"Error: Right fragment {event.right_uid} information not found")
            success = False
    
    # 3. Test get reconstruction results
    _, reconstructed = tree.donors("test")
    
    # Define expected reconstruction results
    recon_expected = {
        "should_have_reconstructions": True  # Expected reconstruction results
    }
    
    if recon_expected["should_have_reconstructions"] and not reconstructed:
        error_messages.append("Error: No reconstruction results obtained")
        success = False
    else:
        print(f"Obtained {len(reconstructed) if reconstructed else 0} reconstruction results")
        # Output reconstruction record information
        if reconstructed:
            for rec in reconstructed:
                recon_type = rec.annotations.get("reconstruction_type", "Unknown")
                original_uid = rec.annotations.get("original_uid", "Unknown")
                seq_len = len(rec.seq)
                print(f"   Reconstruction: Type={recon_type}, OriginalUID={original_uid}, Sequence Length={seq_len}")
    
    # 4. Generate DOT visualization
    dot_str = event_journal.to_graphviz_dot()
    
    # Define expected DOT visualization results
    dot_expected = {
        "contains_nodes": True,  # Expected DOT to contain node information
        "contains_events": True  # Expected DOT to contain event information
    }
    
    # Check if DOT string contains node and event information
    if dot_expected["contains_nodes"] and "node_" not in dot_str:
        error_messages.append("Error: DOT visualization should contain node information")
        success = False
    
    if dot_expected["contains_events"] and "event_" not in dot_str:
        error_messages.append("Error: DOT visualization should contain event information")
        success = False
    
    if success:
        print("✓ Multiple cutting fragment distinction function test passed!")
    else:
        print("✗ Multiple cutting fragment distinction function test failed!")
        for error in error_messages:
            print(f"  {error}")
    
    return success


# If run as main program, execute all tests
if __name__ == "__main__":
    print("Executing tests...")
    
    # Execute multiple cutting test
    success1 = test_multiple_cuts()
    if not success1:
        print("  - Multiple cutting test failed")
        
    # Execute comprehensive nesting test
    success2 = test_comprehensive_nesting()
    if not success2:
        print("  - Comprehensive nesting test failed") 
        
    # Execute multiple cutting fragment distinction function test
    success3 = test_multiple_cuts_fragments_distinction()
    if not success3:
        print("  - Multiple cutting fragment distinction function test failed")
    
    # Display overall test results
    print("\n===== Overall test results =====")
    if success1 and success2 and success3:
        print("✓ All tests passed!")
    else:
        print("✗ Some tests failed!")
        if not success1:
            print("  - Multiple cutting test failed")
        if not success2:
            print("  - Comprehensive nesting test failed") 
        if not success3:
            print("  - Multiple cutting fragment distinction function test failed")
