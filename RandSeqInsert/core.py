#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Tianyu (Sky) Lu (tianyu@lu.fm)

from typing import Tuple, List, Dict, Optional, Set, Union, Iterable, Any
from Bio.SeqRecord import SeqRecord

from utils import create_sequence_record, generate_TSD


class SequenceNode:
    """
    Node in a sequence tree structure, used for efficient sequence insertion operations.
    Implemented as an AVL tree to maintain balance during insertions.
    """
    def __init__(self, data: str, is_donor: bool = False, donor_id: str = None, uid: int = None):
        """
        Initialize a sequence node.

        Args:
            data (str): The sequence string
            is_donor (bool): Whether this node contains a donor sequence
            donor_id (str): Identifier for the donor sequence (if is_donor is True)
            uid (int): Unique identifier for this node
        """
        self.data = data
        self.length = len(data)
        self.is_donor = is_donor
        self.donor_id = donor_id
        self.left = None
        self.right = None
        self.uid = uid

        # Total length of the subtree for efficient traversal
        self.total_length = self.length
        # Height of the node for AVL balancing
        self.height = 1

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
        Convert the tree to a string by in-order traversal.

        Returns:
            str: The concatenated sequence
        """
        return "".join([node.data for node in self])

    def update_height(self):
        """
        Update the height of this node.
        """
        left_height = self.left.height if self.left else 0
        right_height = self.right.height if self.right else 0
        self.height = max(left_height, right_height) + 1

    def update_total_length(self):
        """
        Update the total length of this subtree.

        This method recalculates the total length of the subtree rooted at this node
        by summing the lengths of the left subtree, the current node, and the right subtree.
        """
        left_length = self.left.total_length if self.left else 0
        right_length = self.right.total_length if self.right else 0
        self.total_length = left_length + self.length + right_length

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

        # Create event journal
        from sequenceeventjournal import SequenceEventJournal
        self.event_journal = SequenceEventJournal(self)

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

    def _create_node(self, data: str, is_donor: bool = False, donor_id: str = None, uid: int = None) -> SequenceNode:
        """
        Create a new SequenceNode with a unique UID.

        Args:
            data (str): Sequence data
            is_donor (bool): Whether this node contains a donor sequence
            donor_id (str): Donor ID for tracking and visualization
            uid (int): Preset UID, if provided

        Returns:
            SequenceNode: The newly created node
        """
        uid = uid or self._get_next_uid()
        node = SequenceNode(data, is_donor, donor_id, uid)
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
            donor_seq (str): Donor sequence to insert
            donor_id (str): Identifier for the donor sequence
            tsd_length (int): Length of Target Site Duplication (TSD) to generate
            debug (bool): Enable debug output

        Returns:
            SequenceNode: New root node after insertion
        """
        # Skip insertion if donor sequence is empty
        if not donor_seq:
            return node

        current = node
        parent_stack = []
        path_directions = []  # Record the path direction from root to current node ('left' or 'right')

        # Create donor node UID
        donor_node_uid = self._get_next_uid()

        # Iteratively find insertion position
        while True:
            node_start = current.left.total_length if current.left else 0
            node_end = node_start + current.length

            # Case 1: Position is in the current node's left subtree
            if abs_position <= node_start:
                if current.left:
                    # Record parent node and direction for backtracking
                    parent_stack.append(current)
                    path_directions.append('left')
                    current = current.left
                else:
                    # Create new left child node - use preset UID
                    donor_node = self._create_node(donor_seq, True, donor_id, donor_node_uid)
                    current.left = donor_node
                    
                    
                    current.update()
                    break

            # Case 2: Position is inside the current node
            elif node_start < abs_position < node_end:
                # Calculate relative position
                rel_pos = abs_position - node_start
                left_data = current.data[:rel_pos]
                right_data = current.data[rel_pos:]

                # Handle TSD generation
                tsd_5 = tsd_3 = ""
                if tsd_length > 0:
                    # Extract source TSD sequence from the original sequence
                    source_tsd_seq = right_data[:min(tsd_length, len(right_data))]

                    # Generate TSD sequences (potentially with mutations)
                    tsd_5, tsd_3 = generate_TSD(source_tsd_seq, tsd_length)

                    # Remove source TSD from right_data as it will be duplicated
                    if len(source_tsd_seq) > 0:
                        right_data = right_data[len(source_tsd_seq):]
                    
                    # Add TSD sequences to the donor sequence
                    donor_seq = tsd_5 + donor_seq + tsd_3

                # Save current node information
                old_node_uid = current.uid
                is_current_donor = current.is_donor
                current_donor_id = current.donor_id

                # Create left and right fragment node UIDs
                left_node_uid = self._get_next_uid()
                right_node_uid = self._get_next_uid()

                # Create left node with preset UID
                new_left = self._create_node(left_data, is_current_donor, current_donor_id, left_node_uid)
                if current.left:
                    new_left.left = current.left
                    new_left.update()

                # Create right node with preset UID
                new_right = self._create_node(right_data, is_current_donor, current_donor_id, right_node_uid)
                if current.right:
                    new_right.right = current.right
                    new_right.update()

                # If current node is donor, record insertion event
                if is_current_donor:
                    # Create TSD information (if any)
                    tsd_info = None
                    if tsd_length > 0:
                        tsd_info = {
                            'length': tsd_length,
                            'tsd_5': tsd_5,
                            'tsd_3': tsd_3
                        }
                    
                    # Record insertion event
                    self.event_journal.record_insertion(
                        donor_uid=donor_node_uid,
                        target_uid=old_node_uid,
                        left_uid=left_node_uid,
                        right_uid=right_node_uid,
                        tsd_info=tsd_info
                    )

                # Release old node UID (Do not reuse it)
                # self._release_uid(old_node_uid)

                # Update node information
                current.data = donor_seq
                current.length = len(donor_seq)
                current.is_donor = True
                current.donor_id = donor_id
                current.uid = donor_node_uid
                self.node_dict[donor_node_uid] = current

                # Set new children
                current.left = new_left
                current.right = new_right
                current.update()
                break

            # Case 3: Position is in the current node's right subtree
            elif abs_position >= node_end:
                if current.right:
                    # Record parent node and direction for backtracking
                    parent_stack.append(current)
                    path_directions.append('right')
                    # Adjust absolute position to fit the right subtree's relative position
                    abs_position -= node_end
                    current = current.right
                else:
                    # Create new right child node - use preset UID
                    donor_node = self._create_node(donor_seq, True, donor_id, donor_node_uid)
                    current.right = donor_node
                    
                    
                    current.update()
                    break
            else:
                raise RuntimeError("[ERROR] Should not reach here")

        # Backtrack and update node heights and total lengths while executing balance
        while parent_stack:
            parent = parent_stack.pop()
            direction = path_directions.pop()

            # Update parent node's child reference
            if direction == 'left':
                # Balance current node using helper method
                current = current.balance()
                parent.left = current
            else:  # 'right'
                # Balance current node using helper method
                current = current.balance()
                parent.right = current

            parent.update()

            # Balance parent node
            current = parent.balance()

        return current

    def _insert_recursive(self, node: SequenceNode, abs_position: int, donor_seq: str,
                          donor_id: str = None, tsd_length: int = 0, debug: bool = False) -> SequenceNode:
        """
        Recursively insert a donor sequence at the absolute position in the tree.

        Args:
            node (SequenceNode): Current node
            abs_position (int): Absolute position in the tree to insert at (0-based)
            donor_seq (str): Donor sequence to insert
            donor_id (str): Identifier for the donor sequence
            tsd_length (int): Length of Target Site Duplication (TSD) to generate
            debug (bool): Enable debug output

        Returns:
            SequenceNode: New node after insertion
        """
        # Skip insertion if donor sequence is empty
        if not donor_seq:
            return node

        # Create donor node UID
        donor_node_uid = self._get_next_uid()

        # Calculate positions in tree
        node_start = node.left.total_length if node.left else 0
        node_end = node_start + node.length

        # Case 1: Position is in the current node's left subtree
        if abs_position <= node_start:
            if node.left:
                node.left = self._insert_recursive(node.left, abs_position, donor_seq, donor_id, tsd_length)
            else:
                # Insert as left child - use preset UID
                donor_node = self._create_node(donor_seq, True, donor_id, donor_node_uid)
                node.left = donor_node
                

            node.update()
            return node.balance()

        # Case 2: Position is inside the current node
        if node_start < abs_position < node_end:
            # Split this node
            rel_pos = abs_position - node_start
            left_data = node.data[:rel_pos]
            right_data = node.data[rel_pos:]

            # Handle TSD generation
            tsd_5 = tsd_3 = ""
            if tsd_length > 0:
                # Extract source TSD sequence from the original sequence
                source_tsd_seq = right_data[:min(tsd_length, len(right_data))]

                # Generate TSD sequences (potentially with mutations)
                tsd_5, tsd_3 = generate_TSD(source_tsd_seq, tsd_length)

                # Remove source TSD from right_data as it will be duplicated
                if len(source_tsd_seq) > 0:
                    right_data = right_data[len(source_tsd_seq):]

                # Add TSD sequences to the donor sequence
                donor_seq = tsd_5 + donor_seq + tsd_3

            # Save current node information
            old_node_uid = node.uid
            is_current_donor = node.is_donor
            current_donor_id = node.donor_id

            # Create left and right fragment node UIDs
            left_node_uid = self._get_next_uid()
            right_node_uid = self._get_next_uid()

            # Create new node with preset UID
            new_left = self._create_node(left_data, is_current_donor, current_donor_id, left_node_uid)
            if node.left:
                new_left.left = node.left
                new_left.update()
                # Balance left subtree
                new_left = new_left.balance()

            # Create right node with preset UID
            new_right = self._create_node(right_data, is_current_donor, current_donor_id, right_node_uid)
            if node.right:
                new_right.right = node.right
                new_right.update()
                # Balance right subtree
                new_right = new_right.balance()

            # If current node is donor, record insertion event
            if is_current_donor:
                # Create TSD information (if any)
                tsd_info = None
                if tsd_length > 0:
                    tsd_info = {
                        'length': tsd_length,
                        'tsd_5': tsd_5,
                        'tsd_3': tsd_3
                    }
                
                # Record insertion event
                self.event_journal.record_insertion(
                    donor_uid=donor_node_uid,
                    target_uid=old_node_uid,
                    left_uid=left_node_uid,
                    right_uid=right_node_uid,
                    tsd_info=tsd_info
                )

            # Release old node UID (Do not reuse it)
            # self._release_uid(old_node_uid)

            # Update node information
            node.data = donor_seq
            node.length = len(donor_seq)
            node.is_donor = True
            node.donor_id = donor_id
            node.uid = donor_node_uid
            self.node_dict[donor_node_uid] = node

            # Set new children
            node.left = new_left
            node.right = new_right
            node.update()

            return node.balance()

        # Case 3: Position is in the current node's right subtree
        if abs_position >= node_end:
            if node.right:
                node.right = self._insert_recursive(node.right, abs_position - node_end, donor_seq, donor_id, tsd_length)
            else:
                # Insert as right child - use preset UID
                donor_node = self._create_node(donor_seq, True, donor_id, donor_node_uid)
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
        Collect all donor nodes and reconstruct nested donors.
        Uses in-order traversal to collect only active nodes in the final sequence,
        and uses event journal for efficient reconstruction.

        Args:
            seq_id (str): Original sequence ID

        Returns:
            Tuple[List[SeqRecord], List[SeqRecord]]:
                - Regular donor records (only those without nested insertions)
                - Reconstructed donor records
        """
        # Get active nodes via in-order traversal
        active_nodes = self.collect_active_nodes()
        
        # Collect donor records from active nodes
        donor_records = []
        abs_position_map = self.event_journal._calculate_absolute_positions()
        
        for node in active_nodes:
            if node.is_donor:
                # Get absolute position information
                start_pos = abs_position_map.get(node.uid, 0)
                
                # Convert to 1-based index for output
                start_pos_1based = start_pos + 1
                
                # Get donor sequence and check for TSD
                donor_seq = node.data
                donor_length = node.length
                
                # Find insertion events where this node is the donor
                for event in self.event_journal.events:
                    if event.donor_uid == node.uid and event.tsd_info:
                        # Extract TSD information
                        tsd_info = event.tsd_info
                        tsd_5 = tsd_info.get('tsd_5', '')
                        tsd_3 = tsd_info.get('tsd_3', '')
                        
                        # Remove TSD sequences from donor
                        if tsd_5 and donor_seq.startswith(tsd_5):
                            donor_seq = donor_seq[len(tsd_5):]
                            donor_length -= len(tsd_5)
                            # Adjust position to account for removed 5' TSD
                            start_pos_1based += len(tsd_5)
                        
                        if tsd_3 and donor_seq.endswith(tsd_3):
                            donor_seq = donor_seq[:-len(tsd_3)]
                            donor_length -= len(tsd_3)
                        
                        # We found the event with this donor, no need to continue
                        break
                
                # Skip donors with zero length after TSD removal
                if donor_length <= 0:
                    continue
                    
                # Create donor ID - now reflects sequence without TSD
                donor_id = f"{seq_id}_{start_pos_1based}_{start_pos_1based + donor_length - 1}-+-{donor_length}"
                if node.donor_id:
                    donor_id += f"-{node.donor_id}"
                
                # Create record with TSD-free sequence
                donor_record = create_sequence_record(donor_seq, donor_id)
                donor_record.annotations["uid"] = node.uid
                donor_record.annotations["position"] = start_pos_1based
                donor_record.annotations["length"] = donor_length
                
                # Add to result list
                donor_records.append(donor_record)
        
        # Reconstruct nested sequences using event journal, passing active_nodes
        reconstructed_records = self.event_journal.reconstruct_donors_to_records(seq_id, active_nodes)
        
        # Get the UIDs of donors that have been targets of insertion events (have nested insertions)
        # Pass active_nodes to ensure consistency
        excluded_uids = self.event_journal.get_reconstructed_donor_uids(active_nodes)
        
        # Filter donor_records to only include those that don't have nested insertions
        if excluded_uids:
            donor_records = [record for record in donor_records
                            if record.annotations.get("uid") not in excluded_uids]
        
        return donor_records, reconstructed_records

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
            event_journal: Event journal instance for fragment detection
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

        # Redesigned color assignment logic for clarity:
        # 1. First set base colors (donor = blue, acceptor = green)
        # 2. If it's a cut fragment, change to pink
        # 3. If it has nesting relationship, change to yellow
        # 4. If conditions 2 and 3 are both true, use purple

        # Check if this is a fragment of a cut donor
        if event_journal.is_fragment(node.uid):
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
                         "Length: ", str(node.length), "\\n",
                         cut_half])

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
