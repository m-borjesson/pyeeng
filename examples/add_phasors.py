from transforms import ThreePhase

def main():
	i = ThreePhase()
	i.p=10

	print(i.to_pnz_dict(polar=True))
	print(i.to_dqz_dict(polar=True))

if __name__ == "__main__":
	main()