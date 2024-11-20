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

        # Orbital Input Pane
        self.orbital_frame = None
        self.orbital_entries = []

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

    def create_fragment_section(self):
        self.create_section("Fragment", "Fragments", self.fragment_frame, self.fragment_entries, self.analyze_fragments)

    def create_orbital_section(self):
        self.create_section("Orbital", "Orbitals", self.orbital_frame, self.orbital_entries, self.analyze_orbitals)

    def create_section(self, section_type, section_title, frame, entries, analyze_func):
        if frame is None:
            frame = ttk.Frame(self.root, padding="10")
            frame.grid(row=1 if section_type == "Fragment" else 2, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        for widget in frame.winfo_children():
            widget.destroy()
        entries.clear()

        ttk.Label(frame, text=f"Description of {section_title}").grid(row=0, column=0, columnspan=2, pady=5)
        desc_entry = ttk.Entry(frame, width=40)
        desc_entry.grid(row=1, column=0, columnspan=2, pady=5)
        ttk.Button(
            frame, text=f"Analyze {section_title}", command=lambda: analyze_func(desc_entry.get())
        ).grid(row=2, column=0, columnspan=2, pady=10)

    def analyze_fragments(self, description):
        self.fragment_description = description
        self.analyze_section(description, "Fragments", self.fragment_frame, self.fragment_entries)

    def analyze_orbitals(self, description):
        self.orbital_description = description
        self.analyze_section(description, "Orbitals", self.orbital_frame, self.orbital_entries)

    def analyze_section(self, description, section_title, frame, entries):
        try:
            num_items = self.calculate_items(description)
            if num_items <= 0:
                raise ValueError(f"Number of {section_title.lower()} must be positive.")
        except ValueError:
            messagebox.showerror("Invalid Input", f"Please enter a valid description of {section_title.lower()}.")
            return

        for widget in frame.winfo_children():
            widget.destroy()

        ttk.Label(frame, text=f"Number of {section_title}: {num_items}").grid(row=0, column=0, columnspan=2, pady=5)
        for i in range(num_items):
            row = i + 1
            ttk.Label(frame, text=f"{section_title[:-1]} {i + 1} Type").grid(row=row, column=0, padx=5, pady=5, sticky=tk.W)
            type_entry = ttk.Entry(frame, width=10)
            type_entry.grid(row=row, column=1, padx=5, pady=5)

            ttk.Label(frame, text="Number").grid(row=row, column=2, padx=5, pady=5, sticky=tk.W)
            number_entry = ttk.Entry(frame, width=10)
            number_entry.grid(row=row, column=3, padx=5, pady=5)

            entries.append((type_entry, number_entry))

    def calculate_items(self, description):
        total_items = 0
        groups = description.split()
        for group in groups:
            if '*' in group:
                match = re.match(r"(\d+)\*(\d+)", group)
                if match:
                    total_items += int(match.group(2))
            else:
                total_items += 1
        return total_items

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
        inputs = {key: entry.get() for key, entry in self.entries.items() if entry.get()}

        fragment_inputs = []
        fragment_inputs.append(self.fragment_description)
        for pane in self.fragment_entries:
            fragment_data = {
                "type": pane[0].get(),
                "number": pane[1].get(),
            }
            fragment_inputs.append(fragment_data)
        inputs["fragments"] = fragment_inputs

        orbital_inputs = []
        orbital_inputs.append(self.orbital_description)
        for pane in self.orbital_entries:
            orbital_data = {
                "type": pane[0].get(),
                "number": pane[1].get(),
            }
            orbital_inputs.append(orbital_data)
        inputs["orbitals"] = orbital_inputs

        input_text = "$ctrl\n"
        input_text += " ".join(f"{key}={value}" for key, value in inputs.items() if key not in ["fragments", "orbitals"])
        input_text += "\n$frag\n"
        for frag in fragment_inputs:
            input_text += f"{frag['type']} {frag['number']}\n"
        input_text += "$orb\n"
        for orb in orbital_inputs:
            input_text += f"{orb['type']} {orb['number']}\n"
        input_text += "$end"

        self.show_result(input_text)


# Run the GUI
if __name__ == "__main__":
    root = tk.Tk()
    app = InputGUI(root)
    root.mainloop()
