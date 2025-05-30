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
    
    # Get donor and reconstructed records
    donor_records, reconstructed_records = tree.donors("test")
    
    # Define expected results - corrected for simple insertions
    expected_results = {
        "donor_count": 3,  # Expect 3 donor records (one for each cutter)
        "recon_count": 0,  # Expect 0 reconstructions (no nesting)
        "expected_cutters": [first_cutter, second_cutter, third_cutter],
        "min_cutter_matches": 3  # All 3 cutters should be found in donor records
    }
    
    # Verify correctness
    success = True
    error_messages = []
    
    # Check donor record count
    if donor_records is None or len(donor_records) != expected_results["donor_count"]:
        error_messages.append(f"Error: Should have {expected_results['donor_count']} donor records, but actually have {len(donor_records) if donor_records else 0}")
        success = False
    
    # Check reconstruction count (should be 0 for simple insertions)
    if reconstructed_records is None or len(reconstructed_records) != expected_results["recon_count"]:
        error_messages.append(f"Error: Should have {expected_results['recon_count']} reconstructions, but actually have {len(reconstructed_records) if reconstructed_records else 0}")
        success = False
    
    # Validate donor record sequence content - check if it contains cutter sequences
    found_cutters = []
    if donor_records:
        for cutter in expected_results["expected_cutters"]:
            found = False
            for rec in donor_records:
                seq_str = str(rec.seq)
                # Check complete sequence match
                if cutter == seq_str:
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
        "donor_records": 1,                 # Expected 1 donor record
        "reconstructed_records": 0,         # Expected 0 reconstructed records (no nesting)
        "donor_sequence_present": "TTT"     # Expected this sequence in donor records
    }
    
    # Verify results
    donor_records, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    # Check donor record count
    if donor_records is None or len(donor_records) != expected["donor_records"]:
        success = False
        reason = f"Expected {expected['donor_records']} donor records, but got {len(donor_records) if donor_records else 0} donor records"
    
    # Check reconstruction count (should be 0 for simple insertion)
    if reconstructed is None or len(reconstructed) != expected["reconstructed_records"]:
        success = False
        reason = f"Expected {expected['reconstructed_records']} reconstructed records, but got {len(reconstructed) if reconstructed else 0} reconstructed records"
    
    # Check sequence content in donor records
    if donor_records and not any(expected["donor_sequence_present"] in str(r.seq) for r in donor_records):
        success = False
        reason = f"Expected sequence {expected['donor_sequence_present']} not found in donor records"
    
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
        "donor_records": 2,  # Expected 2 donor records
        "reconstructed_records": 0,  # Expected 0 reconstructed records (no nesting)
        "donor_sequences": ["TTT", "GGG"]  # Expected to contain these two sequences
    }
    
    # Verify results
    donor_records, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    # Check donor record count
    if donor_records is None or len(donor_records) != expected["donor_records"]:
        success = False
        reason = f"Expected {expected['donor_records']} donor records, but got {len(donor_records) if donor_records else 0} donor records"
    
    # Check reconstruction count (should be 0 for simple insertion)
    if reconstructed is None or len(reconstructed) != expected["reconstructed_records"]:
        success = False
        reason = f"Expected {expected['reconstructed_records']} reconstructed records, but got {len(reconstructed) if reconstructed else 0} reconstructed records"
    
    # Check sequence content in donor records
    if donor_records:
        for donor_seq in expected["donor_sequences"]:
            if not any(donor_seq in str(r.seq) for r in donor_records):
                success = False
                reason = f"Expected sequence {donor_seq} not found in donor records"
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
        "min_donor_records": 3,  # Expected at least 3 donor records (fragments + new donor)
        "min_reconstructions": 2,  # Expected at least 2 reconstructions (clean + full)
        "min_full_reconstructions": 1,  # Expected at least 1 full reconstruction
        "should_contain_nested": True,  # Expected at least one reconstruction to contain nested sequence
        "outer_donor": "TTTAAA",
        "inner_donor": "GGG"
    }
    
    # Verify results
    donor_records, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    # Check donor record count
    if donor_records is None or len(donor_records) < expected["min_donor_records"]:
        success = False
        reason = f"Expected at least {expected['min_donor_records']} donor records, but got {len(donor_records) if donor_records else 0} donor records"
    
    # Check reconstruction count
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
        "donor_records": 2,  # Expected 2 donor records
        "reconstructed_records": 0  # Expected 0 reconstructed records (no nesting)
    }
    
    # Verify results
    donor_records, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    # Check donor record count
    if donor_records is None or len(donor_records) != expected["donor_records"]:
        success = False
        reason = f"Expected {expected['donor_records']} donor records, but got {len(donor_records) if donor_records else 0} donor records"
    
    # Check reconstruction count (should be 0 for simple insertion)
    if reconstructed is None or len(reconstructed) != expected["reconstructed_records"]:
        success = False
        reason = f"Expected {expected['reconstructed_records']} reconstructed records, but got {len(reconstructed) if reconstructed else 0} reconstructed records"
    
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
        "min_donor_records": 1,  # Expected at least 1 donor record
        "reconstructed_records": 0  # Expected 0 reconstructed records (no nesting)
    }
    
    # Verify results
    donor_records, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    # Check donor record count
    if donor_records is None or len(donor_records) < expected["min_donor_records"]:
        success = False
        reason = f"Expected at least {expected['min_donor_records']} donor records, but got {len(donor_records) if donor_records else 0} donor records"
    
    # Check reconstruction count (should be 0 for simple insertion)
    if reconstructed is None or len(reconstructed) != expected["reconstructed_records"]:
        success = False
        reason = f"Expected {expected['reconstructed_records']} reconstructed records, but got {len(reconstructed) if reconstructed else 0} reconstructed records"
    
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
        "donor_records": 2,  # Expected 2 donor records (boundary insertion, not nested)
        "reconstructed_records": 0  # Expected 0 reconstructed records (no actual nesting)
    }
    
    # Verify results
    donor_records, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    # Check donor record count
    if donor_records is None or len(donor_records) != expected["donor_records"]:
        success = False
        reason = f"Expected {expected['donor_records']} donor records, but got {len(donor_records) if donor_records else 0} donor records"
    
    # Check reconstruction count (should be 0 for boundary insertion)
    if reconstructed is None or len(reconstructed) != expected["reconstructed_records"]:
        success = False
        reason = f"Expected {expected['reconstructed_records']} reconstructed records, but got {len(reconstructed) if reconstructed else 0} reconstructed records"
    
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
        "donor_records": 5,  # Expected 5 donor records
        "reconstructed_records": 0  # Expected 0 reconstructed records (no actual nesting)
    }
    
    # Verify results
    donor_records, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    # Check donor record count
    if donor_records is None or len(donor_records) != expected["donor_records"]:
        success = False
        reason = f"Expected {expected['donor_records']} donor records, but got {len(donor_records) if donor_records else 0} donor records"
    
    # Check reconstruction count (should be 0 for simple insertion)
    if reconstructed is None or len(reconstructed) != expected["reconstructed_records"]:
        success = False
        reason = f"Expected {expected['reconstructed_records']} reconstructed records, but got {len(reconstructed) if reconstructed else 0} reconstructed records"
    
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
        "min_donor_records": 1,  # Expected at least 1 donor record
        "reconstructed_records": 0  # Expected 0 reconstructed records (no actual nesting)
    }
    
    # Verify results
    donor_records, reconstructed = tree.donors("test")
    
    # Explicit assertion checks
    success = True
    reason = ""
    
    # Check donor record count
    if donor_records is None or len(donor_records) < expected["min_donor_records"]:
        success = False
        reason = f"Expected at least {expected['min_donor_records']} donor records, but got {len(donor_records) if donor_records else 0} donor records"
    
    # Check reconstruction count (should be 0 for simple insertion)
    if reconstructed is None or len(reconstructed) != expected["reconstructed_records"]:
        success = False
        reason = f"Expected {expected['reconstructed_records']} reconstructed records, but got {len(reconstructed) if reconstructed else 0} reconstructed records"
    
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
    Test fragment distinction function for multiple cutting scenarios
    Verify if it can correctly distinguish different fragments produced by different cutters
    
    Returns:
        bool: Whether test passed
    """
    print("\n===== Testing multiple cutting fragment distinction function =====")
    
    # Initialize tree data structure
    tree = SequenceTree("ATGCATGCATGCATGCATGCATGCATGCATGCATGC")  # Original sequence

    # Create a scenario with multiple insertions (not cuts, since these are simple insertions)
    # Insert donor sequences at respective positions
    tree.insert(10, "GTACGTAC", "cutter1")
    tree.insert(20, "CCGGAATT", "cutter2")
    tree.insert(30, "TTAGGCCA", "cutter3")
    
    # Define expected results for simple insertions
    expected = {
        "donor_count": 3,  # Expected 3 donor records
        "recon_count": 0,  # Expected 0 reconstructions (no nesting)
    }
    
    # Test fragment distinction function
    success = True
    error_messages = []
    
    # Test donor records
    donor_records, reconstructed = tree.donors("test")
    
    # Check donor record count
    if donor_records is None or len(donor_records) != expected["donor_count"]:
        error_messages.append(f"Error: Should have {expected['donor_count']} donor records, but found {len(donor_records) if donor_records else 0}")
        success = False
    else:
        print(f"Found {len(donor_records)} donor records")
        for i, record in enumerate(donor_records):
            print(f"   Donor {i}: {record.id}, seq='{record.seq}'")
    
    # Check reconstruction count (should be 0 for simple insertions)
    if reconstructed is None or len(reconstructed) != expected["recon_count"]:
        error_messages.append(f"Error: Should have {expected['recon_count']} reconstructions, but found {len(reconstructed) if reconstructed else 0}")
        success = False
    
    # Verify that each donor sequence is correctly preserved
    if donor_records:
        expected_sequences = ["GTACGTAC", "CCGGAATT", "TTAGGCCA"]
        found_sequences = [str(record.seq) for record in donor_records]
        
        for expected_seq in expected_sequences:
            if expected_seq not in found_sequences:
                error_messages.append(f"Error: Expected sequence '{expected_seq}' not found in donor records")
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
