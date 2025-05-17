#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from core import SequenceTree

# Create a tree with insertions
tree = SequenceTree("ACTGACTGACTG")
print("=== Initial state ===")
print(f"Root UID: {tree.root.uid}")
print(f"node_dict keys: {list(tree.node_dict.keys())}")
print(f"Events: {tree.event_journal.events}")
print("Nodes:")
for uid, node in tree.node_dict.items():
    print(f"  Node UID {uid}: is_donor={node.is_donor}, donor_id={node.donor_id}, seq={node.data[:10]}")

# Insert a donor
print("\n=== After first insertion ===")
tree.insert(5, "DONOR1", "D1")
print(f"Root UID: {tree.root.uid}")
print(f"node_dict keys: {list(tree.node_dict.keys())}")
print(f"Events: {tree.event_journal.events}")
print("Nodes:")
for uid, node in tree.node_dict.items():
    print(f"  Node UID {uid}: is_donor={node.is_donor}, donor_id={node.donor_id}, seq={node.data[:10]}")

# Insert another donor
print("\n=== After second insertion ===")
tree.insert(2, "DONOR2", "D2")
print(f"Root UID: {tree.root.uid}")
print(f"node_dict keys: {list(tree.node_dict.keys())}")
print(f"Events: {tree.event_journal.events}")
print("Nodes:")
for uid, node in tree.node_dict.items():
    print(f"  Node UID {uid}: is_donor={node.is_donor}, donor_id={node.donor_id}, seq={node.data[:10]}")

# Let's see if node 3 is considered a donor
print("\n=== Testing donor status ===")
print(f"node 3 exists: {3 in tree.node_dict}")
if 3 in tree.node_dict:
    print(f"node 3 is_donor flag: {tree.node_dict[3].is_donor}")
    
print(f"journal.is_donor(3): {tree.event_journal.is_donor(3)}")
print(f"Events with donor_uid=3: {[e for e in tree.event_journal.events if e.donor_uid == 3]}") 