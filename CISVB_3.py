import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import re


class InputGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Input File Creator")

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

    def create_input_fields(self):
        ttk.Label(self.frame, text="Enter Input Keywords").grid(
            row=0, column=0, columnspan=2, pady=5
        )
        for idx, key in enumerate(self.keywords):
            ttk.Label(self.frame, text=key).grid(row=idx + 1, column=0, sticky=tk.W, padx=5, pady=5)
            entry = ttk.Entry(self.frame, width=20)
            entry.grid(row=idx + 1, column=1, padx=5, pady=5)
            self.entries[key] = entry

    def create_buttons(self):
        fragment_button = ttk.Button(self.frame, text="Fragments", command=self.create_fragment_section)
        fragment_button.grid(row=len(self.keywords) + 2, column=0, columnspan=2, pady=10)

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

    def create_fragment_section(self):
        # Create a new frame for fragment inputs if it doesn't exist
        if self.fragment_frame is None:
            self.fragment_frame = ttk.Frame(self.root, padding="10")
            self.fragment_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Clear existing fragment panes
        for widget in self.fragment_frame.winfo_children():
            widget.destroy()
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
        self.description=description
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
#                    total_fragments += int(match.group(1)) * int(match.group(2))
                    total_fragments += int(match.group(2))
            else:
                # Handle individual groups
                total_fragments += 1
        return total_fragments

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
        inputs = {key: entry.get() for key, entry in self.entries.items() if entry.get()}

        # Collect fragment inputs if available
        fragment_inputs = []
        fragment_inputs.append(self.description)
        for pane in self.fragment_entries:
            fragment_data = {
                "type": pane[0].get(),
                "number": pane[1].get(),
            }
            fragment_inputs.append(fragment_data)
        inputs["fragments"] = fragment_inputs

        # Generate input string
        input_text = "$ctrl\n"
        input_text += " ".join(f"{key}={value}" for key, value in inputs.items() if key != "fragments")
        input_text += "\n$frag\n"
        for frag in fragment_inputs:
            input_text += f"{frag['type']} {frag['number']}\n"
        input_text += "$end"

        self.show_result(input_text)


# Run the GUI
if __name__ == "__main__":
    root = tk.Tk()
    app = InputGUI(root)
    root.mainloop()

