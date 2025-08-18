from src.int.ascii_art import fox_banner, fox_confused
from src.tools.tilt_angle import run_tilt_angle

def main():
    print(fox_banner)
    print("Select a tool to run:\n")
    print("1) Tilt angle tool")
    print("2) Other tool (coming soon)")
    print("0) Exit")

    try:
        choice = int(input("\nEnter your choice: "))
        if choice == 1:
            run_tilt_angle()
        elif choice == 0:
            print("Goodbye! 🦊")
        else:
            print(fox_confused)
            print("Invalid choice.")
    except Exception as e:
        print(fox_confused)
        print(f"ERROR: {e}")

if __name__ == "__main__":
    main()
|
