#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Tianyu (Sky) Lu (tianyu@lu.fm)

from typing import Tuple, List, Dict, Optional, Set, Union, Iterable, Any
from Bio.SeqRecord import SeqRecord

from utils import create_sequence_record, generate_TSD
from nestingeventjournal import NestingEventJournal


class SequenceNode:
    """
    Node in a sequence tree structure, used for efficient sequence insertion operations.
    Implemented as an AVL tree to maintain balance during insertions.
    """
    def __init__(self, data: str, is_donor: bool = False, donor_id: str = None, uid: int = None,
                 tsd_5: str = "", tsd_3: str = ""):
        """
        Initialize a sequence node.

        Args:
            data (str): The sequence string (without TSD)
            is_donor (bool): Whether this node contains a donor sequence
            donor_id (str): Identifier for the donor sequence (if is_donor is True)
            uid (int): Unique identifier for this node
            tsd_5 (str): 5' TSD sequence
            tsd_3 (str): 3' TSD sequence
        """
        self.data = data
        self.length = len(data)
        self.is_donor = is_donor
        self.donor_id = donor_id
        self.tsd_5 = tsd_5
        self.tsd_3 = tsd_3
        self.left = None
        self.right = None
        self.uid = uid

        # Total length of the subtree for efficient traversal (including TSD)
        self.total_length = self.length + len(tsd_5) + len(tsd_3)
        # Height of the node for AVL balancing
        self.height = 1

    def get_full_sequence(self) -> str:
        """Get the full sequence including TSD."""
        return self.tsd_5 + self.data + self.tsd_3

    def get_clean_sequence(self) -> str:
        """Get the clean sequence without TSD."""
        return self.data

    def has_tsd(self) -> bool:
        """Check if this node has TSD."""
        return bool(self.tsd_5 or self.tsd_3)

    def get_tsd_length(self) -> int:
        """Get total TSD length."""
        return len(self.tsd_5) + len(self.tsd_3)

    def __iter__(self):
        """
        Implement in-order traversal of the tree, iterating all nodes in left-root-right order
        Yields:
            SequenceNode: Each node in the tree
        """
        yield from SequenceNode._inorder_traversal(self)

    @staticmethod
    def _inorder_traversal(node):
        if not node:
            return
        yield from SequenceNode._inorder_traversal(node.left)
        yield node
        yield from SequenceNode._inorder_traversal(node.right)

    def __str__(self) -> str:
        """
        Convert the tree to a string by in-order traversal, including TSD.

        Returns:
            str: The concatenated sequence with TSD
        """
        return "".join([node.get_full_sequence() for node in self])

    def update_height(self):
        """
        Update the height of this node.
        """
        left_height = self.left.height if self.left else 0
        right_height = self.right.height if self.right else 0
        self.height = max(left_height, right_height) + 1

    def update_total_length(self):
        """
        Update the total length of this subtree (including TSD).
        """
        left_length = self.left.total_length if self.left else 0
        right_length = self.right.total_length if self.right else 0
        self.total_length = left_length + self.length + len(self.tsd_5) + len(self.tsd_3) + right_length

    def update(self):
        self.update_height()
        self.update_total_length()

    def get_balance_factor(self):
        """
        Calculate the balance factor of this node.

        Returns:
            int: Balance factor (left height - right height)
        """
        left_height = self.left.height if self.left else 0
        right_height = self.right.height if self.right else 0
        return left_height - right_height

    def rotate_right(self):
        """
        Perform a right rotation on this node.

        Returns:
            SequenceNode: The new root after rotation
        """
        # Store the left child as the new root
        new_root = self.left

        # The left child's right subtree becomes this node's left subtree
        self.left = new_root.right

        # This node becomes the new root's right child
        new_root.right = self

        self.update()
        new_root.update()
        return new_root

    def rotate_left(self):
        """
        Perform a left rotation on this node.

        Returns:
            SequenceNode: The new root after rotation
        """
        # Store the right child as the new root
        new_root = self.right

        # The right child's left subtree becomes this node's right subtree
        self.right = new_root.left

        # This node becomes the new root's left child
        new_root.left = self

        self.update()
        new_root.update()
        return new_root

    def balance(self):
        """
        Balance this node if needed.

        Returns:
            SequenceNode: The new root after balancing
        """
        self.update_height()
        balance = self.get_balance_factor()

        # Left-Left case
        if balance > 1 and (self.left and self.__get_left_balance() >= 0):
            return self.rotate_right()

        # Left-Right case
        if balance > 1 and (self.left and self.__get_left_balance() < 0):
            self.left = self.left.rotate_left()
            return self.rotate_right()

        # Right-Right case
        if balance < -1 and (self.right and self.__get_right_balance() <= 0):
            return self.rotate_left()

        # Right-Left case
        if balance < -1 and (self.right and self.__get_right_balance() > 0):
            self.right = self.right.rotate_right()
            return self.rotate_left()

        return self

    def __get_left_balance(self):
        """Get balance factor of left child"""
        return self.left.get_balance_factor() if self.left else 0

    def __get_right_balance(self):
        """Get balance factor of right child"""
        return self.right.get_balance_factor() if self.right else 0


class SequenceTree:
    """
    Manages a tree of SequenceNode objects, with its own UID management system.
    Provides high-level operations for sequence insertion and traversal.
    """
    def __init__(self, initial_seq: str, base_uid: int = 0):
        """
        Initialize a new sequence tree with a root node containing the initial sequence.

        Args:
            initial_seq (str): The initial sequence to store in the root node
            base_uid (int): Base UID for this tree's UID management system
        """
        # Initialize UID management system
        self.next_uid = base_uid
        self.available_uids = []
        self.node_dict = {}

        # Create the root node
        self.root = self._create_node(initial_seq, False)
        
        # Create event journal for tracking insertion events
        self.event_journal = NestingEventJournal(self)

    def __str__(self) -> str:
        """
        Convert the tree to a string.

        Returns:
            str: The concatenated sequence
        """
        return str(self.root) if self.root else ""

    def __iter__(self):
        yield from self.root

    def _get_next_uid(self, reuse: bool = False) -> int:
        """Get the next available UID from the UID management system"""
        if reuse and self.available_uids:
            return self.available_uids.pop(0)

        uid = self.next_uid
        self.next_uid += 1
        return uid

    def _release_uid(self, uid: int):
        """Release an UID back to the UID management system"""

        # TODO This part should only be used in remove node (future work)
        # if uid in self.node_dict:
        #     del self.node_dict[uid]

        if self.next_uid == uid + 1:
            self.next_uid -= 1
        elif uid not in self.available_uids:
            self.available_uids.append(uid)

    def _create_node(self, data: str, is_donor: bool = False, donor_id: str = None, uid: int = None,
                     tsd_5: str = "", tsd_3: str = "") -> SequenceNode:
        """
        Create a new SequenceNode with a unique UID.

        Args:
            data (str): Sequence data (without TSD)
            is_donor (bool): Whether this node contains a donor sequence
            donor_id (str): Donor ID for tracking and visualization
            uid (int): Preset UID, if provided
            tsd_5 (str): 5' TSD sequence
            tsd_3 (str): 3' TSD sequence

        Returns:
            SequenceNode: The newly created node
        """
        uid = uid or self._get_next_uid()
        node = SequenceNode(data, is_donor, donor_id, uid, tsd_5, tsd_3)
        self.node_dict[uid] = node
        return node

    def insert(self, abs_position: int, donor_seq: str, donor_id: str = None, tsd_length: int = 0,
               recursive: bool = False, debug: bool = False) -> None:
        """
        Insert a donor sequence at the specified position.

        Args:
            abs_position (int): Absolute position for insertion (1-based)
            donor_seq (str): Donor sequence to insert
            donor_id (str): Identifier for the donor sequence
            tsd_length (int): Length of Target Site Duplication (TSD) to generate
            recursive (bool): Whether to use recursive insertion method
            debug (bool): Enable debug output
        """
        # Skip insertion if donor sequence is empty
        if not donor_seq:
            return

        # Convert 1-based position to 0-based for internal processing
        zero_based_position = abs_position - 1

        if not recursive:
            self.root = self._insert_iterative(self.root, zero_based_position, donor_seq, donor_id, tsd_length, debug)
        else:
            self.root = self._insert_recursive(self.root, zero_based_position, donor_seq, donor_id, tsd_length, debug)

    def _insert_iterative(self, node: SequenceNode, abs_position: int, donor_seq: str,
                          donor_id: str = None, tsd_length: int = 0, debug: bool = False) -> SequenceNode:
        """
        Iteratively insert a donor sequence at the absolute position in the tree.

        Args:
            node (SequenceNode): Current root node
            abs_position (int): Absolute position in the tree to insert at (0-based)
            donor_seq (str): Donor sequence to insert (clean, without TSD)
            donor_id (str): Identifier for the donor sequence
            tsd_length (int): Length of Target Site Duplication (TSD) to generate
            debug (bool): Enable debug output

        Returns:
            SequenceNode: New root node after insertion
        """
        if not donor_seq:
            return node

        current = node
        parent_stack = []
        path_directions = []
        donor_node_uid = self._get_next_uid()

        while True:
            node_start = current.left.total_length if current.left else 0
            node_end = node_start + current.length + current.get_tsd_length()

            # Case 1: Left subtree
            if abs_position <= node_start:
                if current.left:
                    parent_stack.append(current)
                    path_directions.append('left')
                    current = current.left
                else:
                    # Boundary insertion at left - extract TSD from current node
                    tsd_5 = tsd_3 = ""
                    if tsd_length > 0:
                        # Extract TSD from beginning of current node's full sequence
                        full_seq = current.get_full_sequence()
                        tsd_source = full_seq[:tsd_length] if tsd_length <= len(full_seq) else full_seq
                        tsd_5, tsd_3 = generate_TSD(tsd_source, len(tsd_source))
                        
                        # Remove TSD from current node
                        tsd_5_len = len(current.tsd_5)
                        if tsd_length <= tsd_5_len:
                            # TSD entirely from 5' TSD
                            current.tsd_5 = current.tsd_5[tsd_length:]
                        elif tsd_length <= tsd_5_len + current.length:
                            # TSD spans 5' TSD and data
                            data_trim = tsd_length - tsd_5_len
                            current.tsd_5 = ""
                            current.data = current.data[data_trim:]
                            current.length = len(current.data)
                        else:
                            # TSD spans entire node content (rare case)
                            tsd_3_trim = tsd_length - tsd_5_len - current.length
                            current.tsd_5 = ""
                            current.data = ""
                            current.length = 0
                            current.tsd_3 = current.tsd_3[tsd_3_trim:] if tsd_3_trim < len(current.tsd_3) else ""
                    
                    # Create donor with TSD
                    donor_node = self._create_node(donor_seq, True, donor_id, donor_node_uid, tsd_5, tsd_3)
                    current.left = donor_node
                    current.update()
                    break

            # Case 2: Inside current node (considering TSD)
            elif node_start < abs_position < node_end:
                # Calculate position within node's full sequence
                rel_pos = abs_position - node_start
                full_seq = current.get_full_sequence()
                
                # Determine where the split occurs
                tsd_5_len = len(current.tsd_5)
                data_len = current.length
                
                # Generate TSD for insertion
                tsd_5 = tsd_3 = ""
                if tsd_length > 0:
                    if rel_pos < len(full_seq):
                        source_seq = full_seq[rel_pos:rel_pos + tsd_length]
                        tsd_5, tsd_3 = generate_TSD(source_seq, tsd_length)

                if rel_pos <= tsd_5_len:
                    # Split in 5' TSD
                    left_tsd5 = current.tsd_5[:rel_pos]
                    right_tsd5 = current.tsd_5[rel_pos:]
                    
                    left_node_uid = self._get_next_uid()
                    right_node_uid = self._get_next_uid()
                    
                    new_left = self._create_node("", current.is_donor, current.donor_id, left_node_uid, left_tsd5, "")
                    new_right = self._create_node(current.data, current.is_donor, current.donor_id, right_node_uid, right_tsd5, current.tsd_3)
                    
                elif rel_pos <= tsd_5_len + data_len:
                    # Split in data
                    data_rel_pos = rel_pos - tsd_5_len
                    left_data = current.data[:data_rel_pos]
                    right_data = current.data[data_rel_pos:]
                    
                    left_node_uid = self._get_next_uid()
                    right_node_uid = self._get_next_uid()
                    
                    new_left = self._create_node(left_data, current.is_donor, current.donor_id, left_node_uid, current.tsd_5, "")
                    new_right = self._create_node(right_data, current.is_donor, current.donor_id, right_node_uid, "", current.tsd_3)
                    
                else:
                    # Split in 3' TSD
                    tsd3_rel_pos = rel_pos - tsd_5_len - data_len
                    left_tsd3 = current.tsd_3[:tsd3_rel_pos]
                    right_tsd3 = current.tsd_3[tsd3_rel_pos:]
                    
                    left_node_uid = self._get_next_uid()
                    right_node_uid = self._get_next_uid()
                    
                    new_left = self._create_node(current.data, current.is_donor, current.donor_id, left_node_uid, current.tsd_5, left_tsd3)
                    new_right = self._create_node("", current.is_donor, current.donor_id, right_node_uid, right_tsd3, "")

                # Set up subtrees
                if current.left:
                    new_left.left = current.left
                    new_left.update()
                if current.right:
                    new_right.right = current.right
                    new_right.update()

                # Save the old node UID and info before transformation
                old_node_uid = current.uid
                is_current_donor = current.is_donor
                current_donor_id = current.donor_id
                original_data = current.data  # Save the data BEFORE any transformation

                # Transform current node to donor with TSD
                current.data = donor_seq
                current.length = len(donor_seq)
                current.is_donor = True
                current.donor_id = donor_id
                current.uid = donor_node_uid
                current.tsd_5 = tsd_5
                current.tsd_3 = tsd_3
                current.left = new_left
                current.right = new_right
                self.node_dict[donor_node_uid] = current
                
                # Record insertion event if inserting into a donor node
                if is_current_donor:
                    # Create TSD info
                    tsd_info = None
                    if tsd_length > 0:
                        tsd_info = {
                            'length': tsd_length,
                            'tsd_5': tsd_5,
                            'tsd_3': tsd_3
                        }
                    
                    # Create a preserved copy of the original target node for journal reference
                    # Use the saved original data, not the transformed data
                    original_target = SequenceNode(original_data, is_current_donor, current_donor_id, old_node_uid)
                    self.node_dict[old_node_uid] = original_target
                    
                    # Record insertion event
                    self.event_journal.record_insertion(
                        donor_uid=donor_node_uid,
                        target_uid=old_node_uid,
                        left_uid=new_left.uid,
                        right_uid=new_right.uid,
                        tsd_info=tsd_info
                    )
                
                current.update()
                break

            # Case 3: Right subtree
            elif abs_position >= node_end:
                if current.right:
                    parent_stack.append(current)
                    path_directions.append('right')
                    abs_position -= node_end
                    current = current.right
                else:
                    # Boundary insertion at right - extract TSD from current node
                    tsd_5 = tsd_3 = ""
                    if tsd_length > 0:
                        # Extract TSD from end of current node's full sequence
                        full_seq = current.get_full_sequence()
                        if tsd_length <= len(full_seq):
                            tsd_source = full_seq[-tsd_length:]
                            start_pos = len(full_seq) - tsd_length
                        else:
                            tsd_source = full_seq
                            start_pos = 0
                        tsd_5, tsd_3 = generate_TSD(tsd_source, len(tsd_source))
                        
                        # Remove TSD from current node
                        tsd_5_len = len(current.tsd_5)
                        data_len = current.length
                        
                        if start_pos >= tsd_5_len + data_len:
                            # TSD entirely from 3' TSD
                            trim_pos = start_pos - tsd_5_len - data_len
                            current.tsd_3 = current.tsd_3[:trim_pos]
                        elif start_pos >= tsd_5_len:
                            # TSD spans data and 3' TSD
                            data_trim_pos = start_pos - tsd_5_len
                            current.data = current.data[:data_trim_pos]
                            current.length = len(current.data)
                            current.tsd_3 = ""
                        else:
                            # TSD spans all parts
                            current.tsd_5 = current.tsd_5[:start_pos]
                            current.data = ""
                            current.length = 0
                            current.tsd_3 = ""
                    
                    # Create donor with TSD
                    donor_node = self._create_node(donor_seq, True, donor_id, donor_node_uid, tsd_5, tsd_3)
                    current.right = donor_node
                    current.update()
                    break

        # Backtrack and balance
        while parent_stack:
            parent = parent_stack.pop()
            direction = path_directions.pop()
            
            if direction == 'left':
                current = current.balance()
                parent.left = current
            else:
                current = current.balance()
                parent.right = current
                
            parent.update()
            current = parent.balance()

        return current

    def _insert_recursive(self, node: SequenceNode, abs_position: int, donor_seq: str,
                          donor_id: str = None, tsd_length: int = 0, debug: bool = False) -> SequenceNode:
        """
        Recursively insert a donor sequence at the absolute position in the tree.

        Args:
            node (SequenceNode): Current node
            abs_position (int): Absolute position in the tree to insert at (0-based)
            donor_seq (str): Donor sequence to insert (clean, without TSD)
            donor_id (str): Identifier for the donor sequence
            tsd_length (int): Length of Target Site Duplication (TSD) to generate
            debug (bool): Enable debug output

        Returns:
            SequenceNode: New node after insertion
        """
        if not donor_seq:
            return node

        donor_node_uid = self._get_next_uid()
        node_start = node.left.total_length if node.left else 0
        node_end = node_start + node.length + node.get_tsd_length()

        # Case 1: Left subtree
        if abs_position <= node_start:
            if node.left:
                node.left = self._insert_recursive(node.left, abs_position, donor_seq, donor_id, tsd_length)
            else:
                # Boundary insertion at left - extract TSD from current node
                tsd_5 = tsd_3 = ""
                if tsd_length > 0:
                    # Extract TSD from beginning of current node's full sequence
                    full_seq = node.get_full_sequence()
                    tsd_source = full_seq[:tsd_length] if tsd_length <= len(full_seq) else full_seq
                    tsd_5, tsd_3 = generate_TSD(tsd_source, len(tsd_source))
                    
                    # Remove TSD from current node
                    tsd_5_len = len(node.tsd_5)
                    if tsd_length <= tsd_5_len:
                        # TSD entirely from 5' TSD
                        node.tsd_5 = node.tsd_5[tsd_length:]
                    elif tsd_length <= tsd_5_len + node.length:
                        # TSD spans 5' TSD and data
                        data_trim = tsd_length - tsd_5_len
                        node.tsd_5 = ""
                        node.data = node.data[data_trim:]
                        node.length = len(node.data)
                    else:
                        # TSD spans entire node content (rare case)
                        tsd_3_trim = tsd_length - tsd_5_len - node.length
                        node.tsd_5 = ""
                        node.data = ""
                        node.length = 0
                        node.tsd_3 = node.tsd_3[tsd_3_trim:] if tsd_3_trim < len(node.tsd_3) else ""
                
                # Create donor with TSD
                donor_node = self._create_node(donor_seq, True, donor_id, donor_node_uid, tsd_5, tsd_3)
                node.left = donor_node
            node.update()
            return node.balance()

        # Case 2: Inside current node
        if node_start < abs_position < node_end:
            rel_pos = abs_position - node_start
            full_seq = node.get_full_sequence()
            tsd_5_len = len(node.tsd_5)
            data_len = node.length

            # Generate TSD
            tsd_5 = tsd_3 = ""
            if tsd_length > 0 and rel_pos < len(full_seq):
                source_seq = full_seq[rel_pos:rel_pos + tsd_length]
                tsd_5, tsd_3 = generate_TSD(source_seq, tsd_length)

            # Split logic
            if rel_pos <= tsd_5_len:
                left_tsd5 = node.tsd_5[:rel_pos]
                right_tsd5 = node.tsd_5[rel_pos:]
                new_left = self._create_node("", node.is_donor, node.donor_id, self._get_next_uid(), left_tsd5, "")
                new_right = self._create_node(node.data, node.is_donor, node.donor_id, self._get_next_uid(), right_tsd5, node.tsd_3)
            elif rel_pos <= tsd_5_len + data_len:
                data_rel_pos = rel_pos - tsd_5_len
                left_data = node.data[:data_rel_pos]
                right_data = node.data[data_rel_pos:]
                new_left = self._create_node(left_data, node.is_donor, node.donor_id, self._get_next_uid(), node.tsd_5, "")
                new_right = self._create_node(right_data, node.is_donor, node.donor_id, self._get_next_uid(), "", node.tsd_3)
            else:
                tsd3_rel_pos = rel_pos - tsd_5_len - data_len
                left_tsd3 = node.tsd_3[:tsd3_rel_pos]
                right_tsd3 = node.tsd_3[tsd3_rel_pos:]
                new_left = self._create_node(node.data, node.is_donor, node.donor_id, self._get_next_uid(), node.tsd_5, left_tsd3)
                new_right = self._create_node("", node.is_donor, node.donor_id, self._get_next_uid(), right_tsd3, "")

            # Set up subtrees
            if node.left:
                new_left.left = node.left
                new_left.update()
                new_left = new_left.balance()
            if node.right:
                new_right.right = node.right
                new_right.update()
                new_right = new_right.balance()

            # Save the old node UID and info before transformation
            old_node_uid = node.uid
            is_current_donor = node.is_donor
            current_donor_id = node.donor_id
            original_data = node.data  # Save the data BEFORE any transformation

            # Transform node to donor
            node.data = donor_seq
            node.length = len(donor_seq)
            node.is_donor = True
            node.donor_id = donor_id
            node.uid = donor_node_uid
            node.tsd_5 = tsd_5
            node.tsd_3 = tsd_3
            node.left = new_left
            node.right = new_right
            self.node_dict[donor_node_uid] = node
            
            # Record insertion event if inserting into a donor node
            if is_current_donor:
                # Create TSD info
                tsd_info = None
                if tsd_length > 0:
                    tsd_info = {
                        'length': tsd_length,
                        'tsd_5': tsd_5,
                        'tsd_3': tsd_3
                    }
                
                # Create a preserved copy of the original target node for journal reference
                # Use the saved original data, not the transformed data
                original_target = SequenceNode(original_data, is_current_donor, current_donor_id, old_node_uid)
                self.node_dict[old_node_uid] = original_target
                
                # Record insertion event
                self.event_journal.record_insertion(
                    donor_uid=donor_node_uid,
                    target_uid=old_node_uid,
                    left_uid=new_left.uid,
                    right_uid=new_right.uid,
                    tsd_info=tsd_info
                )
            
            node.update()
            return node.balance()

        # Case 3: Right subtree
        if abs_position >= node_end:
            if node.right:
                node.right = self._insert_recursive(node.right, abs_position - node_end, donor_seq, donor_id, tsd_length)
            else:
                # Boundary insertion at right - extract TSD from current node
                tsd_5 = tsd_3 = ""
                if tsd_length > 0:
                    # Extract TSD from end of current node's full sequence
                    full_seq = node.get_full_sequence()
                    if tsd_length <= len(full_seq):
                        tsd_source = full_seq[-tsd_length:]
                        start_pos = len(full_seq) - tsd_length
                    else:
                        tsd_source = full_seq
                        start_pos = 0
                    tsd_5, tsd_3 = generate_TSD(tsd_source, len(tsd_source))
                    
                    # Remove TSD from current node
                    tsd_5_len = len(node.tsd_5)
                    data_len = node.length
                    
                    if start_pos >= tsd_5_len + data_len:
                        # TSD entirely from 3' TSD
                        trim_pos = start_pos - tsd_5_len - data_len
                        node.tsd_3 = node.tsd_3[:trim_pos]
                    elif start_pos >= tsd_5_len:
                        # TSD spans data and 3' TSD
                        data_trim_pos = start_pos - tsd_5_len
                        node.data = node.data[:data_trim_pos]
                        node.length = len(node.data)
                        node.tsd_3 = ""
                    else:
                        # TSD spans all parts
                        node.tsd_5 = node.tsd_5[:start_pos]
                        node.data = ""
                        node.length = 0
                        node.tsd_3 = ""
                
                # Create donor with TSD
                donor_node = self._create_node(donor_seq, True, donor_id, donor_node_uid, tsd_5, tsd_3)
                node.right = donor_node
            node.update()
            return node.balance()

        raise RuntimeError("[ERROR] Should not reach here")

    def collect_active_nodes(self) -> List[SequenceNode]:
        """
        Perform an in-order traversal of the tree to collect all active nodes in the final sequence.
        This ensures we only consider nodes that are actually part of the final sequence,
        not those that might have been created but later replaced or discarded.

        Returns:
            List[SequenceNode]: List of all active nodes in the tree, in in-order traversal order
        """
        active_nodes = []
        
        def _inorder_traverse(node):
            if not node:
                return
            _inorder_traverse(node.left)
            active_nodes.append(node)
            _inorder_traverse(node.right)
            
        _inorder_traverse(self.root)
        return active_nodes

    def donors(self, seq_id: str) -> Tuple[List[SeqRecord], List[SeqRecord]]:
        """
        Collect all donor nodes and reconstruct cut donors using tree analysis.

        Args:
            seq_id (str): Original sequence ID

        Returns:
            Tuple[List[SeqRecord], List[SeqRecord]]:
                - Regular donor records 
                - Reconstructed donor records (donors that were cut by others)
        """
        donor_records = []
        reconstructed_records = []
        
        # Collect all donor nodes with position info
        def _collect_donors_with_positions(node, current_pos=0):
            if not node:
                return current_pos
            
            # Process left subtree
            left_end_pos = _collect_donors_with_positions(node.left, current_pos)
            
            # Process current node
            node_start_pos = left_end_pos
            if node.is_donor:
                # Calculate clean sequence position (excluding TSD)
                clean_start = node_start_pos + len(node.tsd_5)
                clean_length = node.length
                clean_end = clean_start + clean_length
                
                # Skip empty donors
                if clean_length > 0:
                    # Create donor ID with 1-based positions
                    donor_id = f"{seq_id}_{clean_start + 1}_{clean_end}-+-{clean_length}"
                    if node.donor_id:
                        donor_id += f"-{node.donor_id}"
                    
                    # Create record with clean sequence
                    donor_record = create_sequence_record(node.get_clean_sequence(), donor_id)
                    donor_record.annotations["uid"] = node.uid
                    donor_record.annotations["position"] = clean_start + 1  # 1-based
                    donor_record.annotations["length"] = clean_length
                    donor_record.annotations["tsd_5"] = node.tsd_5
                    donor_record.annotations["tsd_3"] = node.tsd_3
                    donor_record.annotations["original_donor_id"] = node.donor_id
                    
                    donor_records.append(donor_record)
            
            # Process right subtree
            right_end_pos = _collect_donors_with_positions(node.right, node_start_pos + node.length + node.get_tsd_length())
            
            return right_end_pos
        
        _collect_donors_with_positions(self.root)
        
        # Get active nodes for consistency
        active_nodes = self.collect_active_nodes()
        
        # Get reconstructed donors from event journal
        event_reconstructed_records = self.event_journal.reconstruct_donors_to_records(seq_id, active_nodes)
        
        # Get UIDs of donors that need reconstruction (have nested insertions)
        excluded_uids = self.event_journal.get_reconstructed_donor_uids(active_nodes)
        
        # Filter donor_records to exclude donors that have nested insertions
        if excluded_uids:
            donor_records = [record for record in donor_records
                            if record.annotations.get("uid") not in excluded_uids]
        
        # Combine both reconstruction methods
        # If event journal has reconstructions, use those preferentially
        if event_reconstructed_records:
            reconstructed_records = event_reconstructed_records
        else:
            # Fall back to tree-based reconstruction for legacy compatibility
            # Analyze for reconstruction: find donors that were cut by others
            donor_groups = {}  # Group by original_donor_id
            for record in donor_records:
                orig_id = record.annotations.get("original_donor_id")
                if orig_id:
                    if orig_id not in donor_groups:
                        donor_groups[orig_id] = []
                    donor_groups[orig_id].append(record)
            
            # Reconstruct cut donors
            for orig_id, fragments in donor_groups.items():
                if len(fragments) > 1:  # Multiple fragments = donor was cut
                    # Sort fragments by position
                    fragments.sort(key=lambda x: x.annotations["position"])
                    
                    # Create clean reconstruction (original sequence only)
                    clean_seq = "".join(str(f.seq) for f in fragments)
                    clean_id = f"{seq_id}_reconstructed_{orig_id}_clean"
                    clean_record = create_sequence_record(clean_seq, clean_id)
                    clean_record.annotations["reconstruction_type"] = "clean"
                    clean_record.annotations["original_donor_id"] = orig_id
                    clean_record.annotations["fragment_count"] = len(fragments)
                    
                    # Create full reconstruction (with nested insertions)
                    full_seq = self._get_full_sequence_between_positions(
                        fragments[0].annotations["position"] - 1,  # Convert to 0-based
                        fragments[-1].annotations["position"] + fragments[-1].annotations["length"] - 1
                    )
                    full_id = f"{seq_id}_reconstructed_{orig_id}_full"
                    full_record = create_sequence_record(full_seq, full_id)
                    full_record.annotations["reconstruction_type"] = "full"
                    full_record.annotations["original_donor_id"] = orig_id
                    full_record.annotations["fragment_count"] = len(fragments)
                    
                    reconstructed_records.extend([clean_record, full_record])
        
        return donor_records, reconstructed_records
    
    def _get_full_sequence_between_positions(self, start_pos: int, end_pos: int) -> str:
        """Extract sequence between two positions including all nested content."""
        def _extract_range(node, current_pos=0):
            if not node:
                return "", current_pos
            
            result = ""
            
            # Process left subtree
            left_result, left_end = _extract_range(node.left, current_pos)
            if self._overlaps_range(current_pos, left_end, start_pos, end_pos):
                result += left_result
            
            # Process current node
            node_start = left_end
            node_end = node_start + node.length + node.get_tsd_length()
            
            if self._overlaps_range(node_start, node_end, start_pos, end_pos):
                # Calculate the portion of this node that falls within range
                extract_start = max(0, start_pos - node_start)
                extract_end = min(node.length + node.get_tsd_length(), end_pos - node_start + 1)
                
                full_node_seq = node.get_full_sequence()
                if extract_start < len(full_node_seq) and extract_end > extract_start:
                    result += full_node_seq[extract_start:extract_end]
            
            # Process right subtree
            right_result, right_end = _extract_range(node.right, node_end)
            if self._overlaps_range(node_end, right_end, start_pos, end_pos):
                result += right_result
            
            return result, right_end
        
        full_seq, _ = _extract_range(self.root)
        return full_seq
    
    def _overlaps_range(self, range1_start: int, range1_end: int, range2_start: int, range2_end: int) -> bool:
        """Check if two ranges overlap."""
        return not (range1_end < range2_start or range2_end < range1_start)

    def to_graphviz_dot(self) -> str:
        """
        Generate a Graphviz DOT representation of the tree structure for visualization.

        Returns:
            str: Graphviz DOT format string
        """
        if not self.root:
            return "digraph SequenceTree { }"

        # Initialize the DOT string with graph declaration
        dot_str = ["digraph SequenceTree {",
                   "  bgcolor=\"#FFFFFF\";",
                   "  node [fontcolor=\"#000000\", shape=box, style=filled];",
                   "  edge [fontcolor=\"#000000\", penwidth=2.0];"]

        # Generate nodes and edges through recursive traversal
        nodes, edges = SequenceTree._build_graphviz_dot_nodes_edges(self.root, self.event_journal)

        # Add all nodes and edges to the DOT string
        for node in nodes:
            dot_str.append(f"  {node}")
        for edge in edges:
            dot_str.append(f"  {edge}")

        dot_str.append('}')
        return '\n'.join(dot_str)

    @staticmethod
    def _build_graphviz_dot_nodes_edges(node: SequenceNode, event_journal, 
                                        abs_pos: int = 0, null_leaf: bool = False) -> tuple:
        """
        Recursively build nodes and edges for Graphviz visualization.

        Args:
            node (SequenceNode): Current node
            event_journal: Event journal for the tree
            abs_pos (int): Current absolute position (0-based)

        Returns:
            tuple: (nodes list, edges list)
        """
        if not node:
            return [], []

        nodes = []
        edges = []

        # Calculate positions
        left_length = node.left.total_length if node.left else 0

        # Calculate start and end positions for display
        start_pos = abs_pos + left_length
        end_pos = start_pos + node.length

        # Convert to 1-based positions for display
        start_pos_1based = start_pos + 1
        end_pos_1based = end_pos

        # Use uid directly as node ID
        node_id = f"node_{node.uid}"

        # Determine node type and color
        if node.is_donor:
            node_type = "Donor"
            fill_color = "lightblue"
        else:
            node_type = "Acceptor"
            fill_color = "lightgreen"

        # Process fragment and nesting information
        cut_half = ""

        # Check if this is a fragment of a cut donor
        if event_journal and event_journal.is_fragment(node.uid):
            # Get fragment info (original_uid, is_left, cutter_uid)
            fragment_info = event_journal.get_fragment_info(node.uid)
            if fragment_info:
                original_uid, is_left, cutter_uid = fragment_info
                half_type = "L" if is_left else "R"
                cut_half = f"Cut: {half_type}, by {cutter_uid}"
                fill_color = "lightpink"  # Cut fragments shown in pink

        # Create node label with position information
        label = "".join([node_type, " | ", str(node.uid), "\\n",
                         str(start_pos_1based), "\\l",
                         str(end_pos_1based), "\\l",
                         "Length: ", str(node.length)])
        
        # Add cut information if present
        if cut_half:
            label += "\\n" + cut_half

        # Add the node to the nodes list
        nodes.append(f'{node_id} [label="{label}", fillcolor="{fill_color}"];')

        def create_null_leaf(node_id: str, direction: str) -> tuple:
            invisible_node_id = f"null_{direction}_{node.uid}"
            node = f'{invisible_node_id} [label="NULL_LEAF", style="invis"];'
            edge = f'{node_id} -> {invisible_node_id} [style="invis"];'
            return node, edge

        # Process left child if exists
        if node.left:
            # Left child should start at the same absolute position as its parent
            left_abs_pos = abs_pos
            left_nodes, left_edges = SequenceTree._build_graphviz_dot_nodes_edges(node.left, event_journal, left_abs_pos)
            nodes.extend(left_nodes)
            edges.extend(left_edges)

            # Add edge from this node to left child using uid
            left_id = f"node_{node.left.uid}"
            edges.append(f'{node_id} -> {left_id} [label="L", color="blue"];')
        elif node.right and null_leaf:
            # If no left child but has right child, add invisible node for balance
            null_left_node, null_left_edge = create_null_leaf(node_id, "left")
            nodes.append(null_left_node)
            edges.append(null_left_edge)

        # Process right child if exists
        if node.right:
            # Right child starts at the end position of the current node
            right_abs_pos = abs_pos + left_length + node.length
            right_nodes, right_edges = SequenceTree._build_graphviz_dot_nodes_edges(node.right, event_journal, right_abs_pos)
            nodes.extend(right_nodes)
            edges.extend(right_edges)

            # Add edge from this node to right child using uid
            right_id = f"node_{node.right.uid}"
            edges.append(f'{node_id} -> {right_id} [label="R", color="red"];')
        elif node.left and null_leaf:
            # If no right child but has left child, add invisible node for balance
            null_right_node, null_right_edge = create_null_leaf(node_id, "right")
            nodes.append(null_right_node)
            edges.append(null_right_edge)

        if not node.left and not node.right and null_leaf:
            null_left_node, null_left_edge = create_null_leaf(node_id, "left")
            null_right_node, null_right_edge = create_null_leaf(node_id, "right")
            nodes.extend([null_left_node, null_right_node])
            edges.extend([null_left_edge, null_right_edge])

        return nodes, edges
