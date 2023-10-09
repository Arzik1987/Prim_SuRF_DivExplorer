import numpy as np

class PRIMdens:
    def __init__(self, alpha=0.05):
        self.alpha = alpha
        self.boxes_ = []  # To store intermediate boxes
        self.densities_ = []  # To store densities of intermediate boxes

    def fit(self, X):
        # Initial box is the unit box
        box = np.vstack((np.zeros(X.shape[1]), np.ones(X.shape[1])))
        n_points = X.shape[0]
        
        for iteration in range(100):
            max_density = -np.inf
            best_cut = None
            best_box = None

            for dim in range(X.shape[1]):
                for direction in [0, 1]:  # 0 for lower bound, 1 for upper bound
                    trial_box = box.copy()
                    current_dim_length = trial_box[1, dim] - trial_box[0, dim]
                    
                    # Adjust the boundary in the dimension such that the resulting box's volume is (1-alpha) times the previous box
                    adjustment = current_dim_length * self.alpha
                    if direction == 0:
                        trial_box[0, dim] += adjustment
                    else:
                        trial_box[1, dim] -= adjustment
                    
                    # Count points in the trial box
                    in_box = np.all((X >= trial_box[0]) & (X <= trial_box[1]), axis=1)
                    n_in_box = np.sum(in_box)
                    
                    # Calculate volume of the trial box
                    volume = np.prod(trial_box[1] - trial_box[0])
                    
                    # Calculate density
                    density = n_in_box / volume if volume > 0 else 0
                    
                    if density > max_density:
                        max_density = density
                        best_cut = (dim, direction)
                        best_box = trial_box

            # If we can't find a cut that improves density or the box contains less than 1/alpha points, break
            if best_cut is None or np.sum(np.all((X >= best_box[0]) & (X <= best_box[1]), axis=1)) < 1/self.alpha:
                break

            # Update the box and record it
            box = best_box
            self.boxes_.append(box)
            self.densities_.append(max_density)

        return self

    def get_boxes(self):
        return self.boxes_

    def get_densities(self):
        return self.densities_
