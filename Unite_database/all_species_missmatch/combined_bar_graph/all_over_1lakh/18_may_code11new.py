import matplotlib.pyplot as plt
import numpy as np

# Data extracted from the images
primer_pairs = [
    "ITS86F (F) - ITS4 (R)",
    "ITS3 (F) - ITS4 (R)",
    "ITS1F (F) - ITS2 (R)",
    "ITS5(F) - ITS4(F)",
    "ITS3(F) - LR3(R)",
    "NS7(F) - ITS2(R)",
]

# Data from first image (red bars)
data1 = [10.11, 10.66, 3.66, 3.88, 2.31, 0.27]

# Data from second image (custom blue bars)
data2 = [13.5, 13.1, 5.7, 4.3, 3.0, 0.5]

# Only keep pairs that have data in both images
filtered_pairs = []
filtered_data1 = []
filtered_data2 = []

for pair, d1, d2 in zip(primer_pairs, data1, data2):
    if d1 > 0 and d2 > 0:  # Only include pairs with data in both
        filtered_pairs.append(pair)
        filtered_data1.append(d1)
        filtered_data2.append(d2)

x = np.arange(len(filtered_pairs))
width = 0.3  # Reduced width for less space between bars

fig, ax = plt.subplots(figsize=(8, 5))  # Smaller figure size

# Create bars with red and custom blue colors
rects1 = ax.bar(x - width/2, filtered_data1, width, label='Primer Presence (%)', color='red')
rects2 = ax.bar(x + width/2, filtered_data2, width, label='Species Coverage ≤2 Mismatches (%)', color='#1ECBE1')

# Add labels, title, etc.
ax.set_xlabel('Primer Pair')
ax.set_ylabel('Percentage (%)')
ax.set_title('Comparison of Primer Presence and Species Coverage')
ax.set_xticks(x)
ax.set_xticklabels(filtered_pairs, rotation=45, ha='right')
ax.set_yticks([0, 5, 10, 15, 20, 25])  # Specific y-axis ticks
ax.set_ylim(0, 25)  # Set y-axis limit to 25
ax.legend()

# Add value labels on top of bars
def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f'{height:.1f}%',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)

plt.tight_layout()  # Adjust layout to prevent label cutoff
plt.savefig('combined_primer_comparison.png', dpi=300)
plt.show()
