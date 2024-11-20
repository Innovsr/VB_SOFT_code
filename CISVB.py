import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import json


class InputGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Input File Creator")

        # Main Frame
        self.frame = ttk.Frame(root, padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Input Fields
        self.entries = {}
        self.keywords = [
            "str", "nao", "nae", "nmul" 
        ]
        self.create_input_fields()

        # Buttons
        self.create_buttons()

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
        generate_button = ttk.Button(self.frame, text="Generate Input", command=self.generate_input)
        generate_button.grid(row=len(self.keywords) + 1, column=0, padx=5, pady=10, sticky=tk.E)

        reset_button = ttk.Button(self.frame, text="Reset", command=self.reset_fields)
        reset_button.grid(row=len(self.keywords) + 1, column=1, padx=5, pady=10, sticky=tk.W)

    def generate_input(self):
        # Collect input values into a dictionary
        inputs = {key: entry.get() for key, entry in self.entries.items() if entry.get()}

        # Save to a file in JSON format
        self.save_inputs_to_file(inputs)

        # Display generated input in text format
        input_text = "$ctrl\n" + " ".join(f"{key}={value}" for key, value in inputs.items()) + "\n$end"
        self.show_result(input_text)
        print('inputs',inputs)

    def reset_fields(self):
        for entry in self.entries.values():
            entry.delete(0, tk.END)

    def save_inputs_to_file(self, inputs):
        # Save dictionary to a JSON file
        filename = "inputs.json"
        with open(filename, "w") as file:
            json.dump(inputs, file, indent=4)
        messagebox.showinfo("Input Saved", f"Inputs saved to {filename}")

    def show_result(self, input_text):
        result_window = tk.Toplevel(self.root)
        result_window.title("Generated Input")
        result_window.geometry("400x300")

        text_box = tk.Text(result_window, wrap=tk.WORD)
        text_box.insert(tk.END, input_text)
        text_box.config(state=tk.DISABLED)
        text_box.pack(expand=True, fill=tk.BOTH, padx=10, pady=10)

        ttk.Button(result_window, text="Close", command=result_window.destroy).pack(pady=5)


# Run the GUI
if __name__ == "__main__":
    root = tk.Tk()
    app = InputGUI(root)
    root.mainloop()
