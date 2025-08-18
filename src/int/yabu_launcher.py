import sys
import os

fox_banner = r"""
              /\_/\  
             ( o.o )   YABU Research Tools
              > ^ <  
========================================
"""

fox_confused = r"""
              /\_/\  
             ( o.o )?   YABU is confused...
              > ^ <  
"""

tools = {
    "1": ("Tilt Angle Tool", "tools/tilt_angle_tool.py"),
    # Add more tools here:
    # "2": ("RDF Calculator", "rdf_tool.py"),
    # "3": ("Distance Analyzer", "distance_tool.py"),
}

def main():
    print(fox_banner)
    print("Available tools:\n")
    for k, (name, _) in tools.items():
        print(f" {k}. {name}")
    print(" q. Quit\n")

    choice = input("Select a tool: ").strip()

    if choice.lower() == "q":
        print("Goodbye!")
        sys.exit(0)

    if choice not in tools:
        print(fox_confused)
        print("ERROR: Invalid choice.")
        sys.exit(1)

    tool_name, tool_file = tools[choice]

    if not os.path.exists(tool_file):
        print(fox_confused)
        print(f"ERROR: The script {tool_file} was not found.")
        sys.exit(1)

    print(f"\nLaunching {tool_name}...\n")
    os.system(f"python {tool_file}")

if __name__ == "__main__":
    main()
