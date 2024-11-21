import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import re


class InputGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Input File Creator")

# initialise disctionary to store control data
        self.ctrl_inputs={}
 
        # Main Frame
        self.frame = ttk.Frame(root, padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Input Fields
        self.entries = {}
        self.keywords = ["str", "nao", "nae", "nmul"]
        self.create_input_fields()

        # Buttons
        self.create_buttons()

        # Fragment Input Pane
        self.fragment_frame = None
        self.fragment_entries = []
        self.description = ""

        self.orbital_frame = None
        self.orbital_entries = []
        self.description_orb = ""

    def create_input_fields(self):
        ttk.Label(self.frame, text="Enter Input Keywords").grid(row=0, column=0, columnspan=2, pady=5)
#        ttk.Label(self.frame, text="file_name").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
#        entry = ttk.Emtry(self.frame, width=20)
#        entry.grid(row=1, column=1, padx=5, pady=5)
#        self.entries["file_name"]=entry
        for idx, key in enumerate(self.keywords):
            ttk.Label(self.frame, text=key).grid(row=idx + 1, column=0, sticky=tk.W, padx=5, pady=5)
            entry = ttk.Entry(self.frame, width=20)
            entry.grid(row=idx + 1, column=1, padx=5, pady=5)
            self.entries[key] = entry

    def create_buttons(self):
        fragment_button = ttk.Button(self.frame, text="Fragments", command=self.create_fragment_section)
        fragment_button.grid(row=len(self.keywords) + 2, column=0, columnspan=2, pady=10)

        orbital_button = ttk.Button(self.frame, text="Orbitals", command=self.create_orbital_section)
        orbital_button.grid(row=len(self.keywords) + 3, column=0, columnspan=2, pady=10)

        generate_button = ttk.Button(self.frame, text="Generate Input", command=self.generate_input)
        generate_button.grid(row=len(self.keywords) + 1, column=0, padx=5, pady=10, sticky=tk.E)

        reset_button = ttk.Button(self.frame, text="Reset", command=self.reset_fields)
        reset_button.grid(row=len(self.keywords) + 1, column=1, padx=5, pady=10, sticky=tk.W)

    def reset_fields(self):
        for entry in self.entries.values():
            entry.delete(0, tk.END)
        if self.fragment_entries:
            for pane in self.fragment_entries:
                for field in pane:
                    field.delete(0, tk.END)
        if self.orbital_entries:
            for pane in self.orbital_entries:
                for field in pane:
                    field.delete(0, tk.END)

################################################################################
#### fragment section starts here : first ##
################################################################################
    def create_fragment_section(self):
        # Create a new frame for fragment inputs if it doesn't exist
        if self.fragment_frame is None:
            self.fragment_frame = tk.Toplevel(self.root, padx=10, pady=10)
            self.fragment_frame.title("Fragment inputs")
            self.fragment_frame.geometry("400x300")

#        # Clear existing fragment panes
#        for widget in self.fragment_frame.winfo_children():
#            widget.destroy()
        self.fragment_entries = []

        # Input for description of fragments
        ttk.Label(self.fragment_frame, text="Description of Fragments").grid(
            row=0, column=0, columnspan=2, pady=5
        )
        desc_entry = ttk.Entry(self.fragment_frame, width=40)
        desc_entry.grid(row=1, column=0, columnspan=2, pady=5)
        ttk.Button(
            self.fragment_frame, text="Analyze Description", command=lambda: self.analyze_fragments(desc_entry.get())
        ).grid(row=2, column=0, columnspan=2, pady=10)

    def analyze_fragments(self, description):
        self.description = description
        try:
            # Calculate the total number of fragments
            num_fragments = self.calculate_fragments(description)
            if num_fragments <= 0:
                raise ValueError("Number of fragments must be positive.")
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter a valid description of fragments.")
            return

        # Clear existing fragment panes
        for widget in self.fragment_frame.winfo_children():
            widget.destroy()

        # Create input panes for each fragment
        ttk.Label(self.fragment_frame, text=f"Number of Fragments: {num_fragments}").grid(
            row=0, column=0, columnspan=2, pady=5
        )
        for i in range(num_fragments):
            row = i + 1
            ttk.Label(self.fragment_frame, text=f"Fragment {i + 1} Type").grid(row=row, column=0, padx=5, pady=5, sticky=tk.W)
            type_entry = ttk.Entry(self.fragment_frame, width=10)
            type_entry.grid(row=row, column=1, padx=5, pady=5)

            ttk.Label(self.fragment_frame, text="Number").grid(row=row, column=2, padx=5, pady=5, sticky=tk.W)
            number_entry = ttk.Entry(self.fragment_frame, width=10)
            number_entry.grid(row=row, column=3, padx=5, pady=5)

            self.fragment_entries.append((type_entry, number_entry))

    def calculate_fragments(self, description):
        """
        Parse the description and calculate the total number of fragments.
        """
        total_fragments = 0
        groups = description.split()
        for group in groups:
            if '*' in group:
                # Handle compact notation like "3*4"
                match = re.match(r"(\d+)\*(\d+)", group)
                if match:
                    total_fragments += int(match.group(2))
            else:
                # Handle individual groups
                total_fragments += 1
        return total_fragments

###################################################################################
########## orbital section starts here :
###################################################################################

    def create_orbital_section(self):
        # Create a new frame for orbital inputs if it doesn't exist
        if self.orbital_frame is None:
            self.orbital_frame = tk.Toplevel(self.root, padx=10, pady=10)
            self.orbital_frame.title("orbital inputs")
            self.orbital_frame.geometry("300x300")

#        # Clear existing fragment panes
#        for widget in self.orbital_frame.winfo_children():
#            widget.destroy()
        self.orbital_entries = []

        # Input for description of fragments
        ttk.Label(self.orbital_frame, text="Description of orbitals").grid(
            row=0, column=0, columnspan=2, pady=5
        )
        desc_entry = ttk.Entry(self.orbital_frame, width=40)
        desc_entry.grid(row=1, column=0, columnspan=2, pady=5)
        ttk.Button(
            self.orbital_frame, text="Analyze Description", command=lambda: self.analyze_orbital(desc_entry.get())
        ).grid(row=2, column=0, columnspan=2, pady=10)


    def analyze_orbital(self, description_orb):
        self.description_orb = description_orb
        try:
            # Calculate the total number of orbitals
            num_orbital = self.calculate_orbital(description_orb)
            if num_orbital <= 0:
                raise ValueError("Number of orbitals must be positive.")
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter a valid description of orbitals.")
            return

        # Clear existing fragment panes
        for widget in self.orbital_frame.winfo_children():
            widget.destroy()

        # Create input panes for each fragment
        ttk.Label(self.orbital_frame, text=f"Number of orbitals: {num_orbital}").grid(
            row=0, column=0, columnspan=2, pady=5
        )
        for i in range(num_orbital):
            row = i + 1
            ttk.Label(self.orbital_frame, text="Number").grid(row=row, column=2, padx=5, pady=5, sticky=tk.W)
            number_entry = ttk.Entry(self.orbital_frame, width=10)
            number_entry.grid(row=row, column=3, padx=5, pady=5)

            self.orbital_entries.append(number_entry)


    def calculate_orbital(self, description_orb):
        """
        Parse the description and calculate the total number of fragments.
        """
        total_orbital = 0
        groups = description_orb.split()
        for group in groups:
            if '*' in group:
                # Handle compact notation like "3*4"
                match = re.match(r"(\d+)\*(\d+)", group)
                if match:
                    total_orbital += int(match.group(2))
            else:
                # Handle individual groups
                total_orbital += 1
        return total_orbital

    def show_result(self, input_text):
        result_window = tk.Toplevel(self.root)
        result_window.title("Generated Input")
        result_window.geometry("400x300")

        text_box = tk.Text(result_window, wrap=tk.WORD)
        text_box.insert(tk.END, input_text)
        text_box.config(state=tk.DISABLED)
        text_box.pack(expand=True, fill=tk.BOTH, padx=10, pady=10)

        ttk.Button(result_window, text="Close", command=result_window.destroy).pack(pady=5)

    def generate_input(self):
        # Collect standard inputs
        for key, entry in self.entries.items():
            try:
                value = entry.get()
                if value:
                    self.ctrl_inputs[key]=value
            except AttributeError:
                raise TypeError(f"Expected tk.Entry, but got {type(entry)} for key {key}")
            print('inputs',self.ctrl_inputs)

        # Collect fragment inputs if available
        fragment_inputs = [self.description]
        for pane in self.fragment_entries:
            fragment_type = pane[0].get()
            fragment_number = pane[1].get()
            if fragment_type and fragment_number:
                fragment_inputs.append((fragment_type, fragment_number))
            else:
                raise ValueError(f"Invalid pane in fragment_entries: {pane}")

        orbital_inputs = [self.description_orb]
        for pane in self.orbital_entries:
            orbital_number = pane.get()
            if orbital_number:
                orbital_inputs.append(orbital_number)

        # Generate input string
        input_text = "$ctrl\n"
        input_text += " ".join(f"{key}={value}" for key, value in self.ctrl_inputs.items() if key != "fragments" and key != "orbital")
        input_text += "\n$end"
        input_text += "\n$frag\n"
        input_text += f"{self.description}\n"
        for frag in fragment_inputs[1:]:
            input_text += f"{frag[0]} {frag[1]}\n"
        input_text += "$end"
        input_text += "\n$orb\n"
        input_text += f"{self.description_orb}\n"
        for orb in orbital_inputs[1:]:
            input_text += f"{orb}\n"
        input_text += "$end"

        self.show_result(input_text)

    def get_data(self):
        return self.ctrl_inputs

# Run the GUI
if __name__ == "__main__":
    root = tk.Tk()
#    root1 = tk.Tk()
    app = InputGUI(root)
    root.mainloop()
    ctrl_data=app.get_data()
    print('ctrl data', ctrl_data)
#    root1.mainloop()
